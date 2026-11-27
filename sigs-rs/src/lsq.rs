//! hmf-common sigs.LeastSquaresFit: non-negative least squares (Lawson and Hanson's NNLS, as translated
//! to Java by Alan Kaminsky), step for step so that the contributions equal Java's to the last bit.
//!
//! Solves A x = b for x >= 0, minimising ||A x - b||^2, for an M x N matrix A (the signatures, one column
//! each) and the M counts b.

use anyhow::{bail, Result};

pub struct LeastSquaresFit {
    m: usize,
    n: usize,
    /// A, row-major (a[l * n + j]); replaced by QA while solving
    a: Vec<f64>,
    /// b; replaced by Qb while solving
    b: Vec<f64>,
    contribs: Vec<f64>,
    index: Vec<usize>,
    nsetp: usize,
    /// the squared norm of the residual vector
    pub normsqr: f64,
    w: Vec<f64>,
    zz: Vec<f64>,
    terms: [f64; 2],
    itmax: usize,
}

const FACTOR: f64 = 0.01;

impl LeastSquaresFit {
    pub fn new(rows: usize, cols: usize) -> Result<LeastSquaresFit> {
        if rows == 0 || cols == 0 {
            let msg = format!("invalid rows({rows}) or cols({cols})");
            hmf_common::error!("{msg}");
            bail!("IllegalArgumentException: {msg}");
        }
        Ok(LeastSquaresFit {
            m: rows,
            n: cols,
            a: vec![0.0; rows * cols],
            b: vec![0.0; rows],
            contribs: vec![0.0; cols],
            index: vec![0; cols],
            nsetp: 0,
            normsqr: 0.0,
            w: vec![0.0; cols],
            zz: vec![0.0; rows],
            terms: [0.0; 2],
            itmax: 3 * cols,
        })
    }

    /// initialise(factors, data): the matrix is copied only if it has M rows and the data only if it has M
    /// values (MatrixUtils.copy and VectorUtils.copyVector skip a copy of another length)
    pub fn initialise(&mut self, factors: &[Vec<f64>], data: &[f64]) {
        if factors.len() == self.m {
            let cols = factors.first().map_or(0, |r| r.len());
            for (i, row) in factors.iter().enumerate() {
                for j in 0..cols {
                    self.a[i * self.n + j] = row[j];
                }
            }
        }
        if data.len() == self.m {
            self.b.copy_from_slice(data);
        }
        self.index.iter_mut().for_each(|x| *x = 0);
        self.w.iter_mut().for_each(|x| *x = 0.0);
        self.zz.iter_mut().for_each(|x| *x = 0.0);
        self.terms = [0.0; 2];
    }

    pub fn contribs(&self) -> &[f64] {
        &self.contribs
    }

    #[inline]
    fn at(&self, l: usize, j: usize) -> f64 {
        self.a[l * self.n + j]
    }
    #[inline]
    fn set(&mut self, l: usize, j: usize, v: f64) {
        self.a[l * self.n + j] = v;
    }

    pub fn solve(&mut self) -> Result<()> {
        let (m, n) = (self.m, self.n);
        let mut iter = 0usize;
        for i in 0..n {
            self.contribs[i] = 0.0;
            self.index[i] = i;
        }
        self.nsetp = 0;
        // jj keeps its value between the loops, as in Java
        let mut jj: usize = 0;

        'main: loop {
            if self.nsetp >= n || self.nsetp >= m {
                break 'main;
            }
            // the dual (negative gradient) vector W
            for iz in self.nsetp..n {
                let j = self.index[iz];
                let mut sm = 0.0;
                for l in self.nsetp..m {
                    sm += self.at(l, j) * self.b[l];
                }
                self.w[j] = sm;
            }

            let mut iz;
            let mut j;
            let up;
            'candidate: loop {
                let mut wmax = 0.0;
                let mut izmax: isize = -1;
                for iz2 in self.nsetp..n {
                    let j2 = self.index[iz2];
                    if self.w[j2] > wmax {
                        wmax = self.w[j2];
                        izmax = iz2 as isize;
                    }
                }
                if wmax <= 0.0 {
                    break 'main;
                }
                iz = izmax as usize;
                j = self.index[iz];

                let asave = self.at(self.nsetp, j);
                let up_try = self.construct_householder(self.nsetp, self.nsetp + 1, j)?;
                let mut unorm = 0.0;
                for l in 0..self.nsetp {
                    unorm += sqr(self.at(l, j));
                }
                unorm = unorm.sqrt();
                if (unorm + self.at(self.nsetp, j).abs() * FACTOR) - unorm > 0.0 {
                    self.zz[..m].copy_from_slice(&self.b[..m]);
                    self.apply_householder_vec(self.nsetp, self.nsetp + 1, j, up_try)?;
                    let ztest = self.zz[self.nsetp] / self.at(self.nsetp, j);
                    if ztest > 0.0 {
                        up = up_try;
                        break 'candidate;
                    }
                }
                self.set(self.nsetp, j, asave);
                self.w[j] = 0.0;
            }

            self.b[..m].copy_from_slice(&self.zz[..m]);
            self.index[iz] = self.index[self.nsetp];
            self.index[self.nsetp] = j;
            self.nsetp += 1;

            for jz in self.nsetp..n {
                jj = self.index[jz];
                self.apply_householder_col(self.nsetp - 1, self.nsetp, j, up, jj)?;
            }
            for l in self.nsetp..m {
                self.set(l, j, 0.0);
            }
            self.w[j] = 0.0;

            // solve the triangular system, the solution in zz
            self.solve_triangular(&mut jj);

            'secondary: loop {
                iter += 1;
                if iter > self.itmax {
                    hmf_common::error!("Too many iterations");
                    return Ok(());
                }
                let mut alpha = 2.0;
                for ip in 0..self.nsetp {
                    let l = self.index[ip];
                    if self.zz[ip] <= 0.0 {
                        let t = -self.contribs[l] / (self.zz[ip] - self.contribs[l]);
                        if alpha > t {
                            alpha = t;
                            jj = ip;
                        }
                    }
                }
                if alpha == 2.0 {
                    break 'secondary;
                }
                for ip in 0..self.nsetp {
                    let l = self.index[ip];
                    self.contribs[l] += alpha * (self.zz[ip] - self.contribs[l]);
                }

                let mut i = self.index[jj];
                'tertiary: loop {
                    self.contribs[i] = 0.0;
                    if jj != self.nsetp - 1 {
                        jj += 1;
                        for j2 in jj..self.nsetp {
                            let ii = self.index[j2];
                            self.index[j2 - 1] = ii;
                            let (a1, a2) = (self.at(j2 - 1, ii), self.at(j2, ii));
                            let g = compute_givens(a1, a2, &mut self.terms);
                            self.set(j2 - 1, ii, g);
                            self.set(j2, ii, 0.0);
                            let cc = self.terms[0];
                            let ss = self.terms[1];
                            for l in 0..n {
                                if l != ii {
                                    let temp = self.at(j2 - 1, l);
                                    let next = self.at(j2, l);
                                    self.set(j2 - 1, l, cc * temp + ss * next);
                                    self.set(j2, l, -ss * temp + cc * next);
                                }
                            }
                            let temp = self.b[j2 - 1];
                            let next = self.b[j2];
                            self.b[j2 - 1] = cc * temp + ss * next;
                            self.b[j2] = -ss * temp + cc * next;
                        }
                    }
                    self.nsetp -= 1;
                    self.index[self.nsetp] = i;

                    // the remaining coefficients of set P should be feasible; any that are not (from
                    // rounding) are moved to set Z
                    jj = 0;
                    while jj < self.nsetp {
                        i = self.index[jj];
                        if self.contribs[i] <= 0.0 {
                            continue 'tertiary;
                        }
                        jj += 1;
                    }
                    break 'tertiary;
                }

                self.zz[..m].copy_from_slice(&self.b[..m]);
                self.solve_triangular(&mut jj);
            }

            for ip in 0..self.nsetp {
                let i = self.index[ip];
                self.contribs[i] = self.zz[ip];
            }
        }

        self.normsqr = 0.0;
        for i in self.nsetp..m {
            self.normsqr += sqr(self.b[i]);
        }
        Ok(())
    }

    /// the back substitution of the Java code, with its use of the jj left by the code before it
    fn solve_triangular(&mut self, jj: &mut usize) {
        for l in 0..self.nsetp {
            let mut ip = self.nsetp - l;
            if l != 0 {
                for ii in 0..ip {
                    self.zz[ii] -= self.at(ii, *jj) * self.zz[ip];
                }
            }
            ip -= 1;
            *jj = self.index[ip];
            self.zz[ip] /= self.at(ip, *jj);
        }
    }

    /// constructHouseholderTransform on column `pivotcol` of A
    fn construct_householder(&mut self, ipivot: usize, i1: usize, pivotcol: usize) -> Result<f64> {
        let m = self.m;
        let mut cl = self.at(ipivot, pivotcol).abs();
        for j in i1..m {
            cl = java_max(self.at(j, pivotcol).abs(), cl);
        }
        if cl <= 0.0 {
            bail!("IllegalArgumentException: NonNegativeLeastSquares.constructHouseholderTransform(): Illegal pivot vector");
        }
        let clinv = 1.0 / cl;
        let mut sm = sqr(self.at(ipivot, pivotcol) * clinv);
        for j in i1..m {
            sm += sqr(self.at(j, pivotcol) * clinv);
        }
        cl *= sm.sqrt();
        if self.at(ipivot, pivotcol) > 0.0 {
            cl = -cl;
        }
        let up = self.at(ipivot, pivotcol) - cl;
        self.set(ipivot, pivotcol, cl);
        Ok(up)
    }

    /// applyHouseholderTransform to column `applycol` of A
    fn apply_householder_col(&mut self, ipivot: usize, i1: usize, pivotcol: usize, up: f64, applycol: usize) -> Result<()> {
        let m = self.m;
        let cl = self.at(ipivot, pivotcol).abs();
        if cl <= 0.0 {
            bail!("IllegalArgumentException: NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot vector");
        }
        let mut b = up * self.at(ipivot, pivotcol);
        if b == 0.0 {
            return Ok(());
        } else if b > 0.0 {
            bail!("IllegalArgumentException: NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot element");
        }
        b = 1.0 / b;
        let mut sm = self.at(ipivot, applycol) * up;
        for i in i1..m {
            sm += self.at(i, applycol) * self.at(i, pivotcol);
        }
        if sm != 0.0 {
            sm *= b;
            let v = self.at(ipivot, applycol) + sm * up;
            self.set(ipivot, applycol, v);
            for i in i1..m {
                let v = self.at(i, applycol) + sm * self.at(i, pivotcol);
                self.set(i, applycol, v);
            }
        }
        Ok(())
    }

    /// applyHouseholderTransform to the vector zz
    fn apply_householder_vec(&mut self, ipivot: usize, i1: usize, pivotcol: usize, up: f64) -> Result<()> {
        let m = self.m;
        let cl = self.at(ipivot, pivotcol).abs();
        if cl <= 0.0 {
            bail!("IllegalArgumentException: NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot vector");
        }
        let mut b = up * self.at(ipivot, pivotcol);
        if b == 0.0 {
            return Ok(());
        } else if b > 0.0 {
            bail!("IllegalArgumentException: NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot element");
        }
        b = 1.0 / b;
        let mut sm = self.zz[ipivot] * up;
        for i in i1..m {
            sm += self.zz[i] * self.at(i, pivotcol);
        }
        if sm != 0.0 {
            sm *= b;
            self.zz[ipivot] += sm * up;
            for i in i1..m {
                self.zz[i] += sm * self.at(i, pivotcol);
            }
        }
        Ok(())
    }
}

#[inline]
fn sqr(x: f64) -> f64 {
    x * x
}

/// Math.max(a, b) (NaN if either is NaN)
#[inline]
fn java_max(a: f64, b: f64) -> f64 {
    hmf_common::java_math::java_max(a, b)
}

/// computeGivensRotation: c and s into terms, sqrt(a^2 + b^2) returned
fn compute_givens(a: f64, b: f64, terms: &mut [f64; 2]) -> f64 {
    if a.abs() > b.abs() {
        let xr = b / a;
        let yr = (1.0 + sqr(xr)).sqrt();
        terms[0] = sign(1.0 / yr, a);
        terms[1] = terms[0] * xr;
        a.abs() * yr
    } else if b != 0.0 {
        let xr = a / b;
        let yr = (1.0 + sqr(xr)).sqrt();
        terms[1] = sign(1.0 / yr, b);
        terms[0] = terms[1] * xr;
        b.abs() * yr
    } else {
        terms[0] = 0.0;
        terms[1] = 1.0;
        0.0
    }
}

#[inline]
fn sign(x: f64, y: f64) -> f64 {
    if y >= 0.0 { x } else { -x }
}
