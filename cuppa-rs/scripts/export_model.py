"""Export the CUPPA classifier (a pickled pycuppa CuppaClassifier) into the binary resource cuppa-rs embeds.

Run once, with pycuppa's own environment (scikit-learn 1.3.0, pandas 2.0, numpy < 2):

    python -I export_model.py cuppa_classifier.38.pickle.gz cuppa_model.38.bin.gz

The resource holds every parameter that pycuppa's prediction (python -m cuppa.predict) uses, and nothing that
only training uses. Arrays are stored with their exact values (IEEE bytes, float32 arrays as float32), so the
prediction in Rust sees the same numbers as numpy does.

Format (little-endian), gzip-compressed:
    magic b"CUPPARSMODEL1\n", u32 entry count, then per entry:
    u16 name length, name (UTF-8), u8 kind, u8 ndim, u64 dims[ndim], payload
    kinds: 1 float64, 2 float32, 3 int64, 4 strings (u32 length + UTF-8 each), 5 one UTF-8 text (u32 length +
    UTF-8), 6 an alias of an earlier string list (u32 length + the earlier entry's name)
    Float payloads are byte-shuffled (all first bytes of the elements, then all second bytes, ...), which
    compresses better; the reader undoes it.
"""

import gzip
import hashlib
import struct
import sys

import numpy as np
import pandas as pd
import sklearn
from sklearn.isotonic import IsotonicRegression

from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.components.passthrough import PassthroughRegression

MAGIC = b"CUPPARSMODEL1\n"


class Writer:
    def __init__(self):
        self.entries = []
        self.names = set()
        self.string_lists = {}

    def _add(self, name, kind, dims, payload):
        assert name not in self.names, name
        self.names.add(name)
        name_b = name.encode("utf-8")
        head = struct.pack("<H", len(name_b)) + name_b + struct.pack("<BB", kind, len(dims))
        head += b"".join(struct.pack("<Q", int(d)) for d in dims)
        self.entries.append(head + payload)

    @staticmethod
    def _shuffle(arr):
        raw = np.ascontiguousarray(arr).view(np.uint8).reshape(-1, arr.dtype.itemsize)
        return raw.T.copy().tobytes()

    def floats(self, name, arr, dtype=None):
        arr = np.asarray(arr)
        if dtype is not None:
            assert arr.dtype == dtype, (name, arr.dtype, dtype)
        if arr.dtype == np.float64:
            kind = 1
        elif arr.dtype == np.float32:
            kind = 2
        else:
            raise TypeError(f"{name}: unexpected dtype {arr.dtype}")
        # C order of the given shape; matrices are passed in the layout the reader expects
        arr = np.ascontiguousarray(arr)
        self._add(name, kind, arr.shape, self._shuffle(arr.astype(arr.dtype.newbyteorder("<"))))

    def ints(self, name, arr):
        arr = np.ascontiguousarray(np.asarray(arr), dtype="<i8")
        self._add(name, 3, arr.shape, arr.tobytes())

    def strings(self, name, values):
        values = [str(v) for v in values]
        key = tuple(values)
        if len(values) > 100 and key in self.string_lists:
            # the same list stored before: an alias to it (feature names are repeated in the model)
            target = self.string_lists[key].encode("utf-8")
            self._add(name, 6, (), struct.pack("<I", len(target)) + target)
            return
        self.string_lists.setdefault(key, name)
        payload = b"".join(struct.pack("<I", len(v.encode("utf-8"))) + v.encode("utf-8") for v in values)
        self._add(name, 4, (len(values),), payload)

    def text(self, name, value):
        b = str(value).encode("utf-8")
        self._add(name, 5, (), struct.pack("<I", len(b)) + b)

    def write(self, path):
        body = MAGIC + struct.pack("<I", len(self.entries)) + b"".join(self.entries)
        # mtime 0: the same model gives the same bytes
        with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0, filename="") as fh:
            fh.write(body)


def check(cond, msg):
    if not cond:
        raise SystemExit("export check failed: " + msg)


def export_lr(w, prefix, lr, dtype, with_means=False):
    check(lr.multi_class == "multinomial", prefix + " multinomial")
    check(lr.fit_intercept, prefix + " intercept")
    check(lr.coef_.shape == (len(lr.classes_), len(lr.feature_names_in_)), prefix + " coef shape")
    w.strings(prefix + ".classes", lr.classes_)
    w.strings(prefix + ".features", lr.feature_names_in_)
    w.floats(prefix + ".coef", lr.coef_, dtype)  # classes x features, row-major as in sklearn
    w.floats(prefix + ".intercept", lr.intercept_, dtype)
    if with_means:
        means = lr.X_means
        check(list(means.index) == list(lr.feature_names_in_), prefix + " means index")
        w.floats(prefix + ".x_means", means.values, dtype)


def export_profiles(w, prefix, profiles, dtype):
    # profiles: features x classes; stored class by class (each class's profile contiguous)
    w.strings(prefix + ".profile_features", profiles.index)
    w.strings(prefix + ".profile_classes", profiles.columns)
    w.floats(prefix + ".profiles", profiles.values.T, dtype)
    # the L2 norm of each class profile, as ProfileSimilarityTransformer's cos_sim computes it at every
    # prediction (B = profiles_.transpose(); np.sqrt(np.sum(B**2, axis=1))): it does not depend on the
    # sample, so numpy's own result is stored rather than reproduced (its summation order depends on the
    # array's memory layout, and float32 profiles are summed in float32)
    B = profiles.transpose()
    norms = np.sqrt(np.sum(B ** 2, axis=1))
    w.floats(prefix + ".profile_norms", np.asarray(norms))


def export_calibrator(w, prefix, cal, classes):
    check(not cal.bypass, prefix + " bypass")
    kinds = []
    for c in classes:
        k = cal.calibrators[c]
        if isinstance(k, PassthroughRegression):
            kinds.append("passthrough")
            continue
        check(isinstance(k, IsotonicRegression), prefix + " isotonic")
        check(k.out_of_bounds == "clip", prefix + " clip")
        check(k.increasing_, prefix + " increasing")
        xt, yt = k.X_thresholds_, k.y_thresholds_
        check(xt.dtype == np.float64 and yt.dtype == np.float64, prefix + " threshold dtypes")
        check(len(xt) >= 2, prefix + " thresholds >= 2")
        check(np.all(np.diff(xt) > 0), prefix + " thresholds sorted")
        check(float(k.X_min_) == xt[0] and float(k.X_max_) == xt[-1], prefix + " X_min/X_max")
        kinds.append("isotonic")
        w.floats(f"{prefix}.cal.{c}.x", xt)
        w.floats(f"{prefix}.cal.{c}.y", yt)
    w.strings(prefix + ".cal_kinds", kinds)
    check(set(cal.calibrators.keys()) == set(classes), prefix + " calibrator classes")


def main(pickle_path, out_path):
    check(sklearn.__version__ == "1.3.0", "scikit-learn 1.3.0 required, found " + sklearn.__version__)
    check(np.__version__.startswith("1."), "numpy < 2 required")

    with open(pickle_path, "rb") as fh:
        raw = fh.read()
    digest = hashlib.sha256(raw).hexdigest()

    clf = CuppaClassifier.from_file(pickle_path)
    w = Writer()

    w.text("source.sha256", digest)
    # a cheap identity of the source file for cuppa-rs: its size and last 8 bytes (gzip's CRC-32 and size)
    w.text("source.size", str(len(raw)))
    w.text("source.tail", raw[-8:].hex())
    w.text("source.versions", f"scikit-learn {sklearn.__version__}, numpy {np.__version__}, pandas {pd.__version__}")
    w.strings("classes", clf.classes_)

    req = clf.get_required_features()
    for k in ["gen_pos", "snv96", "event", "sig", "gene_exp", "alt_sj"]:
        w.strings("required." + k, req[k])

    sub = clf["sub_clfs"]
    names = [n for n, _, _ in sub.transformers]
    check(names == ["gen_pos", "snv96", "event", "gene_exp", "alt_sj"], "sub-classifier names")
    for n, _, cols in sub.transformers:
        check(cols.pattern == "^" + n and cols.dtype_include is None and cols.dtype_exclude is None, n + " selector")

    # gen_pos ----------------------------------------------------------------
    gp = clf.gen_pos_clf
    check([s for s, _ in gp.steps] == ["add_noise", "cos_sim", "non_best_scale", "logistic_regression"], "gen_pos steps")
    check(gp["add_noise"].noise_counts == 500, "gen_pos noise counts")
    w.floats("gen_pos.noise_profile", gp["add_noise"].noise_profile, np.float64)
    w.ints("gen_pos.noise_counts", [gp["add_noise"].noise_counts])
    export_profiles(w, "gen_pos", gp["cos_sim"].profiles_, np.float64)
    w.ints("gen_pos.exponent", [gp["non_best_scale"].exponent])
    export_lr(w, "gen_pos.lr", gp["logistic_regression"], np.float64)

    # snv96 ------------------------------------------------------------------
    sn = clf.snv96_clf
    check([s for s, _ in sn.steps] == ["add_noise", "log1p", "logistic_regression"], "snv96 steps")
    w.floats("snv96.noise_profile", sn["add_noise"].noise_profile, np.float64)
    w.ints("snv96.noise_counts", [sn["add_noise"].noise_counts])
    export_lr(w, "snv96.lr", sn["logistic_regression"], np.float64)

    # event ------------------------------------------------------------------
    ev = clf.event_clf
    check([s for s, _ in ev.steps] == ["feature_transform", "logistic_regression"], "event steps")
    ft = ev["feature_transform"]
    check([n for n, _, _ in ft.transformers] == ["other", "sv-tmb"], "event transformers")
    other, svtmb = ft["other"], ft["sv-tmb"]
    check(ft.transformers[0][2].pattern == "^event[.](?:driver|fusion|virus|trait)", "event other pattern")
    check(ft.transformers[1][2].pattern == "^event[.](?:sv|tmb)", "event sv-tmb pattern")
    w.strings("event.other_selected", other.selected_features)
    check([s for s, _ in svtmb.steps] == ["log1p", "max_scale"], "event sv-tmb steps")
    ms = svtmb["max_scale"]
    check(ms.clip, "event max scaler clip")
    w.strings("event.svtmb_features", ms.max_values.index)
    w.floats("event.svtmb_max", ms.max_values.values, np.float64)
    export_lr(w, "event.lr", ev["logistic_regression"], np.float64, with_means=True)

    # RNA sub-classifiers ------------------------------------------------------
    for n, dtype in [("gene_exp", None), ("alt_sj", None)]:
        p = clf[f"sub_clfs"][n]
        check([s for s, _ in p.steps] == ["na_filter", "chi2", "cos_sim", "non_best_scaler", "logistic_regression"], n + " steps")
        na = p["na_filter"]
        check(na.pattern == "^" + n and na.use_first_col, n + " na filter")
        chi2 = p["chi2"]
        check(chi2.clip, n + " max scaler clip")
        check(list(chi2.max_scaler.max_values.index) == list(chi2.selected_features), n + " max index")
        w.strings(n + ".selected", chi2.selected_features)
        mv = chi2.max_scaler.max_values.values
        w.floats(n + ".max", mv)
        prof = p["cos_sim"].profiles_
        export_profiles(w, n, prof, None)
        w.ints(n + ".exponent", [p["non_best_scaler"].exponent])
        export_lr(w, n + ".lr", p["logistic_regression"], None)

    # meta classifiers ---------------------------------------------------------
    meta = clf["meta_clfs"]
    check([n for n, _, _ in meta.transformers] == ["dna_combined", "rna_combined"], "meta names")
    check(meta.transformers[0][2].pattern == "^gen_pos|snv96|event", "dna_combined selector")
    check(meta.transformers[1][2].pattern == "^gene_exp|alt_sj", "rna_combined selector")

    dna = clf.dna_combined_clf
    check([s for s, _ in dna.steps] == ["logistic_regression", "calibrator", "fusion_overrider", "sex_filter"], "dna_combined steps")
    export_lr(w, "dna_combined.lr", dna["logistic_regression"], np.float64)
    export_calibrator(w, "dna_combined", dna["calibrator"], list(dna["logistic_regression"].classes_))
    fo = dna["fusion_overrider"]
    check(not fo.bypass and fo.overrides is not None, "fusion overrider")
    w.floats("dna_combined.mask_base_value", [fo.mask_base_value])
    for col in ["feat_prefix", "feat_basename", "target_class"]:
        w.strings("dna_combined.overrides." + col, fo.overrides[col].astype(str))
    sf = dna["sex_filter"]
    w.strings("dna_combined.male_keywords", sf.keywords_male_classes)
    w.strings("dna_combined.female_keywords", sf.keywords_female_classes)

    rna = clf.rna_combined_clf
    check([s for s, _ in rna.steps] == ["na_filter", "logistic_regression", "calibrator", "sex_filter"], "rna_combined steps")
    check(rna["na_filter"].pattern == "gene_exp|alt_sj" and rna["na_filter"].use_first_col, "rna na filter")
    export_lr(w, "rna_combined.lr", rna["logistic_regression"], np.float64)
    export_calibrator(w, "rna_combined", rna["calibrator"], list(rna["logistic_regression"].classes_))
    sf = rna["sex_filter"]
    w.strings("rna_combined.male_keywords", sf.keywords_male_classes)
    w.strings("rna_combined.female_keywords", sf.keywords_female_classes)

    pc = clf["prob_combiner"]
    check(pc.combine_mode == "multiply", "combiner mode")
    w.floats("combiner.prob_floor", [pc.prob_floor])

    # signature quantiles -------------------------------------------------------
    sq = clf.sig_quantile_transformer
    check(not sq.clip_upper and not sq.clip_lower and sq.include_feat_values, "sig quantile options")
    w.strings("sig.features", sq.feature_names_in_)
    sq_classes = list(sq.transformers_.keys())
    w.strings("sig.classes", sq_classes)
    for c in sq_classes:
        t = sq.transformers_[c]
        check(t.output_distribution == "uniform" and not t.ignore_implicit_zeros, "quantile transformer options")
        check(list(t.feature_names_in_) == list(sq.feature_names_in_), "quantile features")
        w.floats(f"sig.{c}.quantiles", t.quantiles_, np.float64)  # n_quantiles x features
        w.floats(f"sig.{c}.references", t.references_, np.float64)

    # cross-validation performance ------------------------------------------------
    perf = pd.DataFrame(clf.cv_performance)
    check(list(perf.columns) == ["class", "clf_name", "n_total", "n_predicted", "n_correct", "recall", "precision"], "cv columns")
    w.strings("cv.class", perf["class"].astype(str))
    w.strings("cv.clf_name", perf["clf_name"].astype(str))
    for col in ["n_total", "n_predicted", "n_correct"]:
        check(perf[col].dtype == np.int64, "cv int column " + col)
        w.ints("cv." + col, perf[col].values)
    for col in ["recall", "precision"]:
        w.floats("cv." + col, perf[col].values.astype(np.float64), np.float64)

    w.write(out_path)
    print(f"wrote {out_path}: {len(w.entries)} entries, source sha256 {digest}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
