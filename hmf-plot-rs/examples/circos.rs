//! Draws a Circos configuration file: `circos <conf> <out.png> [--ignore-skip-run]`.
fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: circos <conf> <out.png> [--ignore-skip-run]");
        std::process::exit(2);
    }
    let opts = hmf_plot::circos::Options { ignore_skip_run: a.iter().any(|x| x == "--ignore-skip-run") };
    let t = std::time::Instant::now();
    let fig = hmf_plot::circos::draw_conf_file(std::path::Path::new(&a[1]), opts).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        std::process::exit(1)
    });
    let t1 = t.elapsed();
    fig.write_png(&a[2]).unwrap();
    eprintln!("{} items, layout {:?}, total {:?}", fig.items.len(), t1, t.elapsed());
}
