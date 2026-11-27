//! Sets COBALT_BUILD_DATE, the build time (UTC, yyyy-MM-dd HH:mm) the version file holds, as the Maven build
//! of the Java tool puts it in cobalt.version. COBALT_BUILD_DATE in the environment overrides it.

fn main() {
    println!("cargo:rerun-if-env-changed=COBALT_BUILD_DATE");
    println!("cargo:rerun-if-changed=build.rs");
    let date = std::env::var("COBALT_BUILD_DATE").unwrap_or_else(|_| {
        let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs() as i64;
        let (days, rem) = (secs.div_euclid(86400), secs.rem_euclid(86400));
        // civil date from days since 1970-01-01 (H. Hinnant's algorithm)
        let z = days + 719468;
        let era = z.div_euclid(146097);
        let doe = z - era * 146097;
        let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let d = doy - (153 * mp + 2) / 5 + 1;
        let m = if mp < 10 { mp + 3 } else { mp - 9 };
        let y = yoe + era * 400 + if m <= 2 { 1 } else { 0 };
        format!("{y:04}-{m:02}-{d:02} {:02}:{:02}", rem / 3600, rem % 3600 / 60)
    });
    println!("cargo:rustc-env=COBALT_BUILD_DATE={date}");
}
