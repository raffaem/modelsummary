source("helpers.R")

# `metrics_from_gof_map()` narrows the statistics requested from {performance}
# to those an explicit `gof_map` whitelist can actually display.

mk <- function(x) modelsummary:::sanitize_gof_map(x)

# not a whitelist -> unchanged
expect_equal(modelsummary:::metrics_from_gof_map(mk(NULL)), "common")
expect_equal(modelsummary:::metrics_from_gof_map(mk("all")), "common")

# whitelist without R2 -> the expensive metric is not requested
expect_equal(
  modelsummary:::metrics_from_gof_map(mk(c("aic", "bic", "nobs"))),
  c("AIC", "BIC")
)
expect_false("R2" %in% modelsummary:::metrics_from_gof_map(mk(c("rmse", "nobs"))))

# whitelist with R2 -> R2 is requested
expect_true("R2" %in% modelsummary:::metrics_from_gof_map(mk(c("r.squared", "nobs"))))

# adjusted R2 is a distinct metric: some classes (lm) do not return it from a
# bare "R2" request, so asking only for "R2" would silently drop the row
expect_equal(modelsummary:::metrics_from_gof_map(mk("adj.r.squared")), "R2_adj")
expect_equal(
  modelsummary:::metrics_from_gof_map(mk(c("r.squared", "adj.r.squared"))),
  c("R2", "R2_adj")
)

# nothing recognisable -> unchanged, rather than guessing nothing is needed
expect_equal(modelsummary:::metrics_from_gof_map(mk(c("nobs", "vcov.type"))), "common")

# everything wanted -> keep the keyword, which may cover more for some classes
expect_equal(
  modelsummary:::metrics_from_gof_map(
    mk(c("aic", "aicc", "bic", "r.squared", "adj.r.squared", "icc", "rmse"))
  ),
  "common"
)

# --- every statistic `"common"` can produce must be recognised -----------------
# Which R2 flavour `"common"` returns depends on the model class. An unmapped
# flavour would be dropped whenever the whitelist also names something we do
# recognise, e.g. gof_map = c("aic", "r2.tjur") on a logistic glm.
requiet("lme4")
set.seed(1)
n <- 300
dat <- data.frame(
  y = rnbinom(n, mu = 4, size = 2),
  b = rbinom(n, 1, .4),
  x = rnorm(n),
  g = factor(sample(letters[1:8], n, TRUE))
)
mods <- list(
  lm = lm(y ~ x, dat),
  logit = glm(b ~ x, dat, family = binomial),
  poisson = glm(y ~ x, dat, family = poisson),
  lmer = suppressMessages(lme4::lmer(y ~ x + (1 | g), data = dat))
)
for (nm in names(mods)) {
  common <- suppressWarnings(performance::model_performance(
    mods[[nm]], metrics = "common", verbose = FALSE))
  common <- colnames(insight::standardize_names(common, style = "broom"))
  expect_true(all(!is.na(modelsummary:::metrics_for_statistic(common))))
}

# the concrete case: a recognised statistic alongside a class-specific R2
expect_true("R2" %in% modelsummary:::metrics_from_gof_map(mk(c("aic", "r2.tjur"))))
expect_true("R2" %in% modelsummary:::metrics_from_gof_map(mk(c("aic", "r2.nagelkerke"))))
expect_true("R2" %in% modelsummary:::metrics_from_gof_map(mk(c("aic", "r2.conditional"))))

# --- end to end: narrowing must never drop a row the table would have shown ---
mod <- lm(mpg ~ hp + drat, mtcars)

for (gm in list(
  c("nobs", "r.squared"),
  c("nobs", "adj.r.squared"),
  c("nobs", "rmse"),
  c("aic", "bic"),
  c("nobs", "r.squared", "adj.r.squared", "rmse", "aic", "bic")
)) {
  narrowed <- modelsummary(mod, output = "data.frame", gof_map = gm)
  # force the old behaviour for comparison
  wide <- modelsummary(mod, output = "data.frame", gof_map = gm, metrics = "common")
  expect_equivalent(narrowed, wide)
}

# a whitelist naming a statistic {performance} never returns is still served
# from the other backend, and is unaffected by the narrowing
tab <- modelsummary(mod, output = "data.frame", gof_map = c("nobs", "rmse"))
expect_true("Num.Obs." %in% tab$term)
expect_true("RMSE" %in% tab$term)
