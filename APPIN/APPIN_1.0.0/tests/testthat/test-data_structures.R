# test-data_structures.R - Tests structures de données
library(testthat)

# =============================================================================
# TEST: Structure centroids
# =============================================================================
test_that("centroids colonnes requises", {
  c <- data.frame(stain_id="peak1", F2_ppm=1.234, F1_ppm=3.456, Volume=12345.67)
  expect_true(all(c("stain_id","F2_ppm","F1_ppm","Volume") %in% names(c)))
})

test_that("stain_id conventions", {
  is_valid <- function(id) grepl("^(peak|man|manual_box|fused)\\d+$", id)
  expect_true(is_valid("peak1")); expect_true(is_valid("man42")); expect_true(is_valid("fused7"))
  expect_false(is_valid("invalid")); expect_false(is_valid("peak")); expect_false(is_valid(""))
})

# =============================================================================
# TEST: Structure boxes
# =============================================================================
test_that("boxes colonnes requises", {
  b <- data.frame(stain_id="p1", xmin=1, xmax=1.5, ymin=3, ymax=3.5)
  expect_true(all(c("stain_id","xmin","xmax","ymin","ymax") %in% names(b)))
})

test_that("box valide min < max", {
  valid <- function(b) b$xmin < b$xmax && b$ymin < b$ymax
  expect_true(valid(data.frame(xmin=1, xmax=2, ymin=3, ymax=4)))
  expect_false(valid(data.frame(xmin=2, xmax=1, ymin=3, ymax=4)))
})

# =============================================================================
# TEST: Pending changes
# =============================================================================
test_that("pending status add/delete", {
  p <- data.frame(stain_id=c("m1","p5"), F2_ppm=c(1,2), F1_ppm=c(3,4), Volume=c(100,200), status=c("add","delete"))
  expect_true(all(p$status %in% c("add","delete")))
})

test_that("apply_pending fonctionne", {
  apply_p <- function(cur, pen) {
    if (is.null(pen) || nrow(pen)==0) return(cur)
    to_add <- pen[pen$status=="add",]; to_del <- pen[pen$status=="delete",]
    if (nrow(to_add) > 0) { to_add$status <- NULL; cur <- rbind(cur, to_add) }
    if (nrow(to_del) > 0) cur <- cur[!cur$stain_id %in% to_del$stain_id,]
    cur
  }
  cur <- data.frame(stain_id=c("p1","p2"), F2_ppm=c(1,2), F1_ppm=c(3,4), Volume=c(100,200))
  pen <- data.frame(stain_id=c("m1","p1"), F2_ppm=c(5,1), F1_ppm=c(6,3), Volume=c(300,100), status=c("add","delete"))
  res <- apply_p(cur, pen)
  expect_equal(nrow(res), 2)
  expect_true("m1" %in% res$stain_id); expect_false("p1" %in% res$stain_id)
})

# =============================================================================
# TEST: Fusion de pics
# =============================================================================
test_that("fusion centroïde pondéré", {
  fuse <- function(p) {
    tot <- sum(p$Volume)
    data.frame(stain_id="fused1", F2_ppm=sum(p$F2_ppm*p$Volume)/tot, F1_ppm=sum(p$F1_ppm*p$Volume)/tot, Volume=tot)
  }
  p <- data.frame(stain_id=c("p1","p2"), F2_ppm=c(1,2), F1_ppm=c(3,4), Volume=c(100,100))
  f <- fuse(p)
  expect_equal(f$F2_ppm, 1.5); expect_equal(f$F1_ppm, 3.5); expect_equal(f$Volume, 200)
})

test_that("fusion poids inégaux", {
  fuse <- function(p) {
    tot <- sum(p$Volume)
    data.frame(stain_id="fused1", F2_ppm=sum(p$F2_ppm*p$Volume)/tot, F1_ppm=sum(p$F1_ppm*p$Volume)/tot, Volume=tot)
  }
  p <- data.frame(stain_id=c("p1","p2"), F2_ppm=c(1,2), F1_ppm=c(3,4), Volume=c(300,100))
  f <- fuse(p)
  expect_equal(f$F2_ppm, 1.25); expect_equal(f$F1_ppm, 3.25)
})

# =============================================================================
# TEST: Génération IDs
# =============================================================================
test_that("generate_new_id unique", {
  gen <- function(pfx, ex) {
    nums <- as.integer(sub(pfx, "", ex[grepl(paste0("^",pfx), ex)]))
    paste0(pfx, if(length(nums)==0) 1 else max(nums)+1)
  }
  ex <- c("peak1","peak2","peak5","man1")
  expect_equal(gen("peak", ex), "peak6")
  expect_equal(gen("man", ex), "man2")
  expect_equal(gen("fused", ex), "fused1")
})

# =============================================================================
# TEST: Paramètres par défaut
# =============================================================================
test_that("Defaults TOCSY", {
  d <- list(contour_start=80000, intensity_threshold=4000, eps=0.0068)
  expect_equal(d$contour_start, 80000)
})

test_that("Defaults HSQC", {
  d <- list(contour_start=20000, intensity_threshold=200, eps=0.002)
  expect_equal(d$contour_start, 20000)
})

test_that("Defaults COSY", {
  d <- list(contour_start=60000, intensity_threshold=20000, eps=0.014)
  expect_equal(d$eps, 0.014)
})

test_that("Defaults UFCOSY", {
  d <- list(contour_start=50000, intensity_threshold=20000, eps=0.014)
  expect_equal(d$contour_start, 50000)
})
