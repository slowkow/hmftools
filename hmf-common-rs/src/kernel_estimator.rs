//! hmf-common utils.KernelEstimator (from Weka, GPL): a weighted kernel density estimate with a normal
//! kernel, values rounded to a precision, as PURPLE's somatic peak search uses it.

/// KernelEstimator(precision, bandwidth)
pub struct KernelEstimator {
    values: Vec<f64>,
    weights: Vec<f64>,
    sum_of_weights: f64,
    standard_dev: f64,
    precision: f64,
}

const MAX_ERROR: f64 = 0.01;

impl KernelEstimator {
    pub fn new(precision: f64, bandwidth: f64) -> KernelEstimator {
        KernelEstimator { values: Vec::new(), weights: Vec::new(), sum_of_weights: 0.0, standard_dev: bandwidth, precision }
    }

    fn round(&self, data: f64) -> f64 {
        (data / self.precision).round_ties_even() * self.precision
    }

    fn find_nearest_value(&self, key: f64) -> usize {
        let (mut low, mut high) = (0usize, self.values.len());
        while low < high {
            let middle = (low + high) / 2;
            let current = self.values[middle];
            if current == key {
                return middle;
            }
            if current > key {
                high = middle;
            } else if current < key {
                low = middle + 1;
            }
        }
        low
    }

    pub fn add_value(&mut self, data: f64, weight: f64) {
        if weight == 0.0 {
            return;
        }
        let data = self.round(data);
        let i = self.find_nearest_value(data);
        if self.values.len() > i && self.values[i] == data {
            self.weights[i] += weight;
        } else {
            self.values.insert(i, data);
            self.weights.insert(i, weight);
        }
        self.sum_of_weights += weight;
    }

    pub fn probability(&self, data: f64) -> f64 {
        let half = self.precision / 2.0;
        if self.values.is_empty() {
            let z_lower = (data - half) / self.standard_dev;
            let z_upper = (data + half) / self.standard_dev;
            return normal_probability(z_upper) - normal_probability(z_lower);
        }
        let mut sum = 0.0;
        let mut weight_sum = 0.0;
        let start = self.find_nearest_value(data);
        for i in start..self.values.len() {
            let delta = self.values[i] - data;
            let p = normal_probability((delta + half) / self.standard_dev) - normal_probability((delta - half) / self.standard_dev);
            sum += p * self.weights[i];
            weight_sum += self.weights[i];
            if p * (self.sum_of_weights - weight_sum) < sum * MAX_ERROR {
                break;
            }
        }
        let mut i = start as i64 - 1;
        while i >= 0 {
            let iu = i as usize;
            let delta = self.values[iu] - data;
            let p = normal_probability((delta + half) / self.standard_dev) - normal_probability((delta - half) / self.standard_dev);
            sum += p * self.weights[iu];
            weight_sum += self.weights[iu];
            if p * (self.sum_of_weights - weight_sum) < sum * MAX_ERROR {
                break;
            }
            i -= 1;
        }
        sum / self.sum_of_weights
    }
}

fn polevl(x: f64, coef: &[f64], n: usize) -> f64 {
    let mut ans = coef[0];
    for c in coef.iter().take(n + 1).skip(1) {
        ans = ans * x + c;
    }
    ans
}

fn p1evl(x: f64, coef: &[f64], n: usize) -> f64 {
    let mut ans = x + coef[0];
    for c in coef.iter().take(n).skip(1) {
        ans = ans * x + c;
    }
    ans
}

pub fn normal_probability(a: f64) -> f64 {
    let x = a * 0.7071067811865476;
    let z = x.abs();
    if z < 0.7071067811865476 {
        0.5 + 0.5 * error_function(x)
    } else {
        let y = 0.5 * error_function_complemented(z);
        if x > 0.0 { 1.0 - y } else { y }
    }
}

pub fn error_function(x: f64) -> f64 {
    const T: [f64; 5] = [9.604973739870516, 90.02601972038427, 2232.005345946843, 7003.325141128051, 55592.30130103949];
    const U: [f64; 5] = [33.56171416475031, 521.3579497801527, 4594.323829709801, 22629.000061389095, 49267.39426086359];
    if x.abs() > 1.0 {
        return 1.0 - error_function_complemented(x);
    }
    let z = x * x;
    x * polevl(z, &T, 4) / p1evl(z, &U, 5)
}

pub fn error_function_complemented(a: f64) -> f64 {
    const P: [f64; 9] = [
        2.461969814735305E-10, 0.5641895648310689, 7.463210564422699, 48.63719709856814, 196.5208329560771, 526.4451949954773, 934.5285271719576,
        1027.5518868951572, 557.5353353693994,
    ];
    const Q: [f64; 8] = [
        13.228195115474499, 86.70721408859897, 354.9377788878199, 975.7085017432055, 1823.9091668790973, 2246.3376081871097, 1656.6630919416134,
        557.5353408177277,
    ];
    const R: [f64; 6] = [0.5641895835477551, 1.275366707599781, 5.019050422511805, 6.160210979930536, 7.4097426995044895, 2.9788666537210022];
    const S: [f64; 6] = [2.2605286322011726, 9.396035249380015, 12.048953980809666, 17.08144507475659, 9.608968090632859, 3.369076451000815];
    let x = a.abs();
    if x < 1.0 {
        return 1.0 - error_function(a);
    }
    let mut z = -a * a;
    if z < -709.782712893384 {
        return if a < 0.0 { 2.0 } else { 0.0 };
    }
    z = crate::java_strictmath::exp(z);
    let (p, q) = if x < 8.0 { (polevl(x, &P, 8), p1evl(x, &Q, 8)) } else { (polevl(x, &R, 5), p1evl(x, &S, 6)) };
    let mut y = z * p / q;
    if a < 0.0 {
        y = 2.0 - y;
    }
    if y == 0.0 {
        if a < 0.0 { 2.0 } else { 0.0 }
    } else {
        y
    }
}
