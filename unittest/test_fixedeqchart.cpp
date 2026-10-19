// Unit tests of model/fixedeqchart.cpp against the oracle's fixtures (docs/agent/CODE_PLAN.md 3.2).
// The pass lines are decision 026's: every entry of Q within a relative max(1e-12, 1e-14 kappa),
// and residuals at most 1e-10. Everything else is compared exactly.

#include "fixedeqchart.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fq = fixedeq;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using fq::Status;

namespace {

const double RESIDUAL = 1e-10;                               // decision 026
const double LO = std::log(1e-5), HI = std::log(10.0);      // the default box (decisions 007, 027)
const double NaN = std::numeric_limits<double>::quiet_NaN();

double qLimit(double kappa) { return std::max(1e-12, 1e-14 * kappa); }

std::string g17(double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

// the oracle's rel_err: the largest per-entry relative difference, the diagonal included
double relErr(const MatrixXd &A, const MatrixXd &B) {
    return ((A - B).array().abs() / B.array().abs()).maxCoeff();
}

double number(const std::string &token) {
    char *end = nullptr;
    double v = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0')
        ADD_FAILURE() << "unreadable number in a fixture: " << token;
    return v;
}

struct Fixture {
    std::map<std::string, std::string> header;
    std::map<std::string, std::vector<double>> arrays;

    bool has(const std::string &label) const { return arrays.count(label) > 0; }
    int n() const { return std::stoi(header.at("n")); }
    VectorXd vec(const std::string &label) const {
        const std::vector<double> &v = arrays.at(label);
        return Eigen::Map<const VectorXd>(v.data(), v.size());
    }
    MatrixXd mat(const std::string &label) const {  // stored row-major
        const std::vector<double> &v = arrays.at(label);
        const int k = n();
        MatrixXd M(k, k);
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
                M(i, j) = v.at(i * k + j);
        return M;
    }
    std::vector<int> ref() const {
        std::vector<int> r;
        for (double v : arrays.at("ref"))
            r.push_back(static_cast<int>(v));
        return r;
    }
};

std::string fixturePath(const std::string &name) { return std::string(FIXTURE_DIR) + "/" + name + ".txt"; }

Fixture readFixture(const std::string &name) {
    std::ifstream in(fixturePath(name));
    if (!in)
        ADD_FAILURE() << "cannot open " << fixturePath(name);
    Fixture f;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("# ", 0) == 0) {
            const size_t colon = line.find(": ", 2);
            f.header[line.substr(2, colon - 2)] = colon == std::string::npos ? "" : line.substr(colon + 2);
        } else if (!line.empty()) {
            std::istringstream s(line);
            std::string label, token;
            s >> label;
            std::vector<double> &values = f.arrays[label];
            while (s >> token)
                values.push_back(number(token));
        }
    }
    return f;
}

std::vector<std::string> fixtureNames() {
    std::vector<std::string> names;
    std::error_code ec;
    for (const auto &entry : std::filesystem::directory_iterator(FIXTURE_DIR, ec))
        if (entry.path().extension() == ".txt" && entry.path().filename() != "targets.txt")
            names.push_back(entry.path().stem().string());
    std::sort(names.begin(), names.end());
    return names;
}

}  // namespace

// ---------------------------------------------------------------- fixtures from the oracle

TEST(Fixtures, AllArePresent) { EXPECT_EQ(fixtureNames().size(), 27u); }

class ChartFixture : public ::testing::TestWithParam<std::string> {};

TEST_P(ChartFixture, MatchesOracle) {
    const Fixture f = readFixture(GetParam());
    const int n = f.n();
    const VectorXd pi = f.vec("pi"), z = f.vec("z");
    const MatrixXd Q_fix = f.mat("Q");
    const std::vector<int> ref = f.ref();
    ASSERT_EQ(pi.size(), n);
    ASSERT_EQ(z.size(), n * (n - 2));
    ASSERT_EQ(static_cast<int>(ref.size()), n);
    const double kappa = number(f.header.at("kappa")), limit = qLimit(kappa);

    // build from the fixture's coordinates
    MatrixXd Q;
    VectorXd nu;
    ASSERT_EQ(fq::build(z, pi, ref, Q, nu), Status::ok);
    const double build_error = relErr(Q, Q_fix);
    EXPECT_LE(build_error, limit);

    // the inverse, checked through build (no separate line on z)
    VectorXd z2;
    MatrixXd Q2;
    VectorXd nu2;
    ASSERT_EQ(fq::inverse(Q_fix, ref, z2), Status::ok);
    ASSERT_EQ(fq::build(z2, pi, ref, Q2, nu2), Status::ok);
    const double roundtrip_error = relErr(Q2, Q_fix);
    EXPECT_LE(roundtrip_error, limit);

    // residuals of both built generators
    const fq::Diagnostics d1 = fq::diagnostics(Q, pi, nu), d2 = fq::diagnostics(Q2, pi, nu2);
    const double stationarity = std::max(d1.stationarity, d2.stationarity);
    const double row_sums = std::max(d1.row_sums, d2.row_sums);
    const double mean_rate_error = std::max(d1.mean_rate_error, d2.mean_rate_error);
    EXPECT_LE(stationarity, RESIDUAL);
    EXPECT_LE(row_sums, RESIDUAL);
    EXPECT_LE(mean_rate_error, RESIDUAL);
    EXPECT_GT(std::min(d1.min_rate, d2.min_rate), 0.0);

    // references
    const std::string rule = f.header.at("ref_rule");
    if (rule == "row_maxima")
        EXPECT_EQ(fq::references(Q_fix), ref);
    if (rule == "lg") {
        MatrixXd lg;
        ASSERT_EQ(fq::reversibleSeed(f.mat("R_ref"), pi, lg), Status::ok);
        EXPECT_EQ(fq::references(lg), ref);
    }

    // seeds
    std::string seed_error = "";
    if (f.has("R")) {
        MatrixXd S;
        ASSERT_EQ(fq::reversibleSeed(f.mat("R"), pi, S), Status::ok);
        EXPECT_LE(relErr(S, Q_fix), limit);
        seed_error = g17(relErr(S, Q_fix));
    }
    if (f.has("Q0")) {
        MatrixXd S;
        ASSERT_EQ(fq::transferSeed(f.mat("Q0"), pi, S), Status::ok);
        EXPECT_LE(relErr(S, Q_fix), limit);
        seed_error = g17(relErr(S, Q_fix));
    }

    // one row of the residual table; googletest 1.12.1 records numbers only as int, hence strings
    RecordProperty("fixture", GetParam());
    RecordProperty("truth", f.header.at("truth"));
    RecordProperty("n", n);
    RecordProperty("kappa", g17(kappa));
    RecordProperty("q_limit", g17(limit));
    RecordProperty("build_error", g17(build_error));
    RecordProperty("roundtrip_error", g17(roundtrip_error));
    RecordProperty("seed_error", seed_error);
    RecordProperty("stationarity", g17(stationarity));
    RecordProperty("row_sums", g17(row_sums));
    RecordProperty("mean_rate_error", g17(mean_rate_error));
    RecordProperty("nu_residual", g17(std::max(d1.nu_residual, d2.nu_residual)));
    RecordProperty("exit_rate_mismatch", g17(std::max(d1.exit_rate_mismatch, d2.exit_rate_mismatch)));
    RecordProperty("min_rate", g17(std::min(d1.min_rate, d2.min_rate)));
}

INSTANTIATE_TEST_SUITE_P(Oracle, ChartFixture, ::testing::ValuesIn(fixtureNames()),
                         [](const ::testing::TestParamInfo<std::string> &info) { return info.param; });

// ---------------------------------------------------------------- references

TEST(References, TiesGoToTheLowerIndex) {
    MatrixXd Q(4, 4);
    Q << -3.0, 1.0, 1.0, 1.0,
          2.0, -4.0, 1.0, 1.0,
          1.0, 3.0, -7.0, 3.0,
          0.5, 0.5, 2.0, -3.0;
    EXPECT_EQ(fq::references(Q), (std::vector<int>{1, 0, 1, 2}));
}

TEST(References, LgAtTheTurtleTargetMatchesTheManifest) {
    const Fixture f = readFixture("lg_turtle");
    MatrixXd Q;
    ASSERT_EQ(fq::reversibleSeed(f.mat("R"), f.vec("pi"), Q), Status::ok);
    // test_scripts/fixedeq/manifest/g0_manifest.json, references.indices (SKDEAEDAYVIRLYATSFFI)
    const std::vector<int> manifest = {15, 11, 3, 6, 0, 6, 3, 0, 18, 19, 9, 1, 10, 18, 0, 16, 15, 13, 13, 9};
    EXPECT_EQ(fq::references(Q), manifest);
}

// ---------------------------------------------------------------- rejections leave outputs untouched

namespace {

struct Outputs {
    MatrixXd Q = MatrixXd::Constant(4, 4, 7.0);
    VectorXd nu = VectorXd::Constant(4, 7.0);
    bool untouched() const { return Q.rows() == 4 && (Q.array() == 7.0).all() && (nu.array() == 7.0).all(); }
};

Status build4(const VectorXd &z, const VectorXd &pi, Outputs &out) {
    return fq::build(z, pi, {1, 0, 0, 0}, out.Q, out.nu);
}

}  // namespace

TEST(Rejections, NonFiniteCoordinate) {
    Outputs out;
    VectorXd z = VectorXd::Zero(8);
    z(0) = NaN;
    EXPECT_EQ(build4(z, VectorXd::Constant(4, 0.25), out), Status::nonFiniteCoordinate);
    EXPECT_TRUE(out.untouched());
}

TEST(Rejections, SoftmaxUnderflowAfterExp) {
    Outputs out;
    VectorXd z = VectorXd::Zero(8);
    z(0) = 1e4;  // exp(0 - 1e4) of the reference destination is exactly 0
    EXPECT_EQ(build4(z, VectorXd::Constant(4, 0.25), out), Status::softmaxUnderflow);
    EXPECT_TRUE(out.untouched());
}

TEST(Rejections, SoftmaxUnderflowAfterDivision) {
    // row 0 (reference 1, coordinates for destinations 2 and 3): exp(-744.5) is the smallest
    // subnormal, nonzero, but divided by the row sum 2 it rounds to 0. The oracle misses this case.
    ASSERT_EQ(std::exp(-744.5), std::numeric_limits<double>::denorm_min());
    Outputs out;
    VectorXd z = VectorXd::Zero(8);
    z(1) = -744.5;
    EXPECT_EQ(build4(z, VectorXd::Constant(4, 0.25), out), Status::softmaxUnderflow);
    EXPECT_TRUE(out.untouched());
}

TEST(Rejections, RateUnderflow) {
    // nu is about (1/3, 4/9, 2/9), so q_02 = (nu_0 / 0.9) * denorm_min rounds to 0
    const double d = std::numeric_limits<double>::denorm_min();
    MatrixXd K(3, 3);
    K << 0.0, 1.0 - d, d,
         0.5, 0.0, 0.5,
         0.5, 0.5, 0.0;
    VectorXd pi(3);
    pi << 0.9, 0.05, 0.05;
    MatrixXd Q = MatrixXd::Constant(3, 3, 7.0);
    VectorXd nu = VectorXd::Constant(3, 7.0);
    EXPECT_EQ(fq::generatorFromJump(K, pi, Q, nu), Status::nonPositiveRate);
    EXPECT_TRUE((Q.array() == 7.0).all() && (nu.array() == 7.0).all());
}

TEST(Rejections, InverseOfANonPositiveOrNonFiniteRate) {
    MatrixXd Q(3, 3);
    Q << -1.0, 0.0, 1.0,
          0.5, -1.0, 0.5,
          0.5, 0.5, -1.0;
    VectorXd z = VectorXd::Constant(3, 7.0);
    EXPECT_EQ(fq::inverse(Q, {2, 0, 0}, z), Status::nonPositiveRate);
    Q(0, 1) = NaN;
    EXPECT_EQ(fq::inverse(Q, {2, 0, 0}, z), Status::nonFiniteRate);
    EXPECT_TRUE((z.array() == 7.0).all());
}

TEST(Rejections, InvalidTarget) {
    Outputs out;
    VectorXd pi(4);
    pi << 0.5, 0.5, 0.0, 0.0;
    EXPECT_EQ(build4(VectorXd::Zero(8), pi, out), Status::invalidTarget);
    pi << 0.5, 0.5, NaN, 0.0;
    EXPECT_EQ(build4(VectorXd::Zero(8), pi, out), Status::invalidTarget);
    EXPECT_TRUE(out.untouched());
}

TEST(Rejections, InvalidInput) {
    Outputs out;
    const VectorXd pi = VectorXd::Constant(4, 0.25);
    EXPECT_EQ(build4(VectorXd::Zero(7), pi, out), Status::invalidInput);
    EXPECT_EQ(fq::build(VectorXd::Zero(8), pi, {0, 0, 0, 0}, out.Q, out.nu), Status::invalidInput);
    EXPECT_EQ(fq::build(VectorXd::Zero(8), pi, {1, 0, 0, 4}, out.Q, out.nu), Status::invalidInput);
    EXPECT_EQ(fq::build(VectorXd::Zero(8), VectorXd::Constant(3, 0.25), {1, 0, 0, 0}, out.Q, out.nu),
              Status::invalidInput);
    EXPECT_TRUE(out.untouched());
}

// ---------------------------------------------------------------- target parser

namespace {

struct TargetCase {
    std::string name, text, verdict;
    double min_freq = 0.0, raw_sum = 0.0;
    std::vector<double> values;
};

std::vector<TargetCase> readTargets() {
    std::ifstream in(fixturePath("targets"));
    std::vector<TargetCase> cases;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("case ", 0) == 0) {
            cases.emplace_back();
            cases.back().name = line.substr(5);
        } else if (line.rfind("text ", 0) == 0) {
            cases.back().text = line.substr(5);
        } else if (line.rfind("min_freq ", 0) == 0) {
            cases.back().min_freq = number(line.substr(9));
        } else if (line.rfind("verdict ", 0) == 0) {
            cases.back().verdict = line.substr(8);
        } else if (line.rfind("raw_sum ", 0) == 0) {
            cases.back().raw_sum = number(line.substr(8));
        } else if (line.rfind("values ", 0) == 0) {
            std::istringstream s(line.substr(7));
            std::string token;
            while (s >> token)
                cases.back().values.push_back(number(token));
        }
    }
    return cases;
}

std::string uniform(int k, const std::string &sep = ",") {
    std::string s;
    for (int i = 0; i < k; i++)
        s += (i ? sep : "") + std::string("0.05");
    return s;
}

}  // namespace

TEST(Target, MatchesTheOracleVerdicts) {
    const std::vector<TargetCase> cases = readTargets();
    ASSERT_EQ(cases.size(), 25u);
    for (const TargetCase &c : cases) {
        SCOPED_TRACE(c.name);
        VectorXd pi = VectorXd::Constant(20, 7.0);
        double raw = 7.0;
        std::string error;
        const Status s = fq::parseTarget(c.text, 20, c.min_freq, pi, raw, error);
        if (c.verdict == "accept") {
            ASSERT_EQ(s, Status::ok) << error;
            ASSERT_EQ(c.values.size(), 20u);
            EXPECT_EQ(raw, c.raw_sum);
            for (int i = 0; i < 20; i++)
                EXPECT_EQ(pi(i), c.values[i]) << "entry " << i;
        } else {
            EXPECT_EQ(s, Status::invalidTarget);
            EXPECT_FALSE(error.empty());
            EXPECT_EQ(raw, 7.0);
            EXPECT_TRUE((pi.array() == 7.0).all());
        }
    }
}

TEST(Target, RefusesNonAsciiInputThatPythonAccepts) {
    // Python's float() reads "0.0_5" and Arabic-Indic digits, and the oracle's \s splits on a
    // no-break space; the compiled parser takes ASCII decimals and ASCII separators only
    const std::string rest = "," + uniform(19);
    for (const std::string &text : {std::string("0.0_5") + rest,
                                    std::string("\xd9\xa0.\xd9\xa0\xd9\xa5") + rest,
                                    std::string("0.05\xc2\xa0") + uniform(19)}) {
        VectorXd pi;
        double raw = 0.0;
        std::string error;
        EXPECT_EQ(fq::parseTarget(text, 20, 1e-4, pi, raw, error), Status::invalidTarget) << text;
    }
}

TEST(Target, CountAndMinimumAreParameters) {
    VectorXd pi;
    double raw = 0.0;
    std::string error;
    ASSERT_EQ(fq::parseTarget("0.5, 0.25, 0.25", 3, 1e-4, pi, raw, error), Status::ok) << error;
    EXPECT_EQ(raw, 1.0);
    EXPECT_EQ(pi, (VectorXd(3) << 0.5, 0.25, 0.25).finished());
    EXPECT_EQ(fq::parseTarget("0.5, 0.25, 0.25", 20, 1e-4, pi, raw, error), Status::invalidTarget);
    EXPECT_EQ(fq::parseTarget("0.5, 0.25, 0.25", 3, 0.3, pi, raw, error), Status::invalidTarget);
}

// ---------------------------------------------------------------- derivative step

namespace {

// Knuth's TwoSum: a + b = s + err exactly in real arithmetic (round to nearest, no contraction)
::testing::AssertionResult exactDisplacement(double x, double h, double probe) {
    const double s = x + h, bv = s - x, av = s - bv;
    const double err = (x - av) + (h - bv);
    if (s == probe && err == 0.0)
        return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure() << "x + h = " << g17(s) << " + " << g17(err) << ", probe " << g17(probe);
}

}  // namespace

// When |x| >= h0, (x + h0) - x is computed without rounding (Sterbenz), so the returned h is
// exactly the distance the probe moved; it can still differ from h0, which rounds when x + h0 is
// formed.
TEST(Step, IsTheProbesDisplacementAndExactWhenXDominates) {
    // no upper bound in the way, so every step is forward
    const double inf = std::numeric_limits<double>::infinity();
    for (double x : {0.0, 1e-8, -1e-8, 1e-3, -1e-3, 0.5, 5.0, -5.0, -11.5, 2.0}) {
        SCOPED_TRACE(x);
        const double h0 = fq::STEP_ETA * std::max(fq::STEP_SCALE, std::fabs(x));
        const double h = fq::scaledStep(x, inf);
        const double probe = x + h0;
        EXPECT_EQ(h, probe - x);
        if (std::fabs(x) >= h0)
            EXPECT_TRUE(exactDisplacement(x, h, probe));
        EXPECT_GT(h, 0.0);
        EXPECT_GE(h, h0 / 2);  // the legacy rule's 1e-4 |x| would give 1e-12 at x = 1e-8
    }
    // the check would catch a step that returned h0 itself: at x = 2 the probe moves by h0 - 2.2e-17
    EXPECT_FALSE(exactDisplacement(2.0, 2e-4, 2.0 + 2e-4));
}

TEST(Step, StepsBackAtTheUpperBound) {
    for (double x : {HI, HI - 1e-5, HI + 1.0}) {
        SCOPED_TRACE(x);
        const double h0 = fq::STEP_ETA * std::max(fq::STEP_SCALE, std::fabs(x));
        const double h = fq::scaledStep(x, HI);
        EXPECT_LT(h, 0.0);
        EXPECT_EQ(h, (x - h0) - x);
        EXPECT_TRUE(exactDisplacement(x, h, x - h0));
    }
    EXPECT_GT(fq::scaledStep(HI - 1e-3, HI), 0.0);
}

TEST(Step, StepsForwardBelowTheLowerBound) {
    // a start outside the box is evaluated as given (decision 027)
    const double x = LO - 5.0, h = fq::scaledStep(x, HI);
    EXPECT_GT(h, 0.0);
    EXPECT_TRUE(exactDisplacement(x, h, x + fq::STEP_ETA * std::fabs(x)));
}

// ---------------------------------------------------------------- diagnostics

TEST(Diagnostics, ReportWhatTheyMeasure) {
    // pi Q = 0 for this Q, mean rate 3/2, nu = pi * exit = (1, 1, 1)/2 unnormalized; every number
    // below is a dyadic rational, so each value is exact
    VectorXd pi(3), nu(3);
    pi << 0.5, 0.25, 0.25;
    nu << 0.5, 0.5, 0.5;
    MatrixXd Q(3, 3);
    Q << -1.0, 0.5, 0.5,
          1.0, -2.0, 1.0,
          1.0, 1.0, -2.0;
    fq::Diagnostics d = fq::diagnostics(Q, pi, nu);
    EXPECT_EQ(d.stationarity, 0.0);
    EXPECT_EQ(d.row_sums, 0.0);
    EXPECT_EQ(d.mean_rate_error, 0.5);
    EXPECT_EQ(d.min_rate, 0.5);
    EXPECT_EQ(d.min_flux, 0.25);
    EXPECT_EQ(d.nu_residual, 0.0);
    EXPECT_EQ(d.exit_rate_mismatch, 0.0);

    Q(0, 1) += 0.125;  // the diagonal is left alone
    d = fq::diagnostics(Q, pi, nu);
    EXPECT_EQ(d.stationarity, 0.0625);   // (pi Q)_1 = pi_0 delta
    EXPECT_EQ(d.row_sums, 0.125);        // row 0 sums to delta; its exit rate is 1
    EXPECT_EQ(d.mean_rate_error, 0.5);   // the diagonal is unchanged
    EXPECT_EQ(d.min_rate, 0.5);
    EXPECT_EQ(d.min_flux, 0.25);
    EXPECT_EQ(d.nu_residual, 0.0625);    // (nu K)_1 - nu_1 = (1/2)(5/8) + (1/2)(1/2) - 1/2
    EXPECT_EQ(d.exit_rate_mismatch, 0.0);
}
