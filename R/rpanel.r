#' @keywords internal
plott.smooth1 <- function(panel) {
   if (panel$method == "manual") panel$h <- panel$h.manual
   panel$opt$se   <- panel$se
   panel$opt$test <- panel$test
   if (panel$model == "none") panel$opt$test <- FALSE
   with(panel, {
      result <- sm.regression.1d(x, y, h, design.mat, 
                  model, weights, rawdata, opt)
      df.h    <- stats::approx(hvec, dfvec, h, rule = 2)$y
      df.text <- paste("df =", as.character(round(df.h, 1)))
      h.text  <- paste("   h =", as.character(signif(h, 3)))
      if ("p" %in% names(result))
         p.text <- paste("   p =", round(result$p, 3))
      else
         p.text <- ""
      graphics::title(paste(df.text, h.text, p.text))
      })
   panel
   }
   
set.bandwidth <- function(panel) {
   if (panel$method != "manual") {
      panel$h <- h.select(panel$x, panel$y, panel$weights, method = panel$method)
      panel$h <- panel$h[1]
      }
   if (is.matrix(panel$x)) ndim <- 2 else ndim <- 1
   if (panel$opt$panel.plot) {
      if (ndim == 1) rpanel::rp.do(panel, replott.smooth1)
      else           rpanel::rp.do(panel, replott.smooth2)
   }
   else {
      if (ndim == 1) rpanel::rp.do(panel, plott.smooth1)
      else           rpanel::rp.do(panel, plott.smooth2)
      }
   
   panel
   }

replott.smooth1 <- function(panel) {
   rpanel::rp.tkrreplot(panel, plot)
   panel
   }

rp.smooth1 <- function(x, y, h, design.mat, model, weights, rawdata, opt) {

   opt <- sm.options(opt)
   replace.na(opt, se,         FALSE)
   replace.na(opt, test,       FALSE)
   replace.na(opt, panel.plot, TRUE)
   opt$verbose <- 0
   nvec  <- 20
   hvec  <- rep(0, nvec)
   hvec[1]    <- h.select(x, y, weights = weights, method = "df", df = 2.1, nbins = 0)
   hvec[nvec] <- h.select(x, y, weights = weights, method = "df", df = 20,  nbins = 0)
   hvec       <- exp(seq(log(hvec[1]), log(hvec[nvec]), length = nvec))
   dfvec      <- seq(2.1, 20, length = nvec)
   for (i in 2:(nvec - 1))
      dfvec[i] <- sum(diag(sm.weight(x, x, hvec[i], weights = weights)))
   if (opt$panel.plot & !requireNamespace("tkrplot", quietly = TRUE)) opt$panel.plot <- FALSE

   smooth.panel <- rpanel::rp.control("Nonparametric regression - 1 covariate", 
                      x = x, y = y, h = h, design.mat = design.mat,
                      model = model, weights = weights, rawdata = rawdata,
                      opt = opt, hvec = hvec, dfvec = dfvec, h.manual = h,
                      method = "manual", se = opt$se, test = opt$test)
   if (opt$panel.plot) {
      rpanel::rp.tkrplot(smooth.panel, plot, plott.smooth1, pos = "right",
                 hscale = opt$hscale, vscale = opt$vscale)
      plotfun <- replott.smooth1
      }
   else
      plotfun <- plott.smooth1
   rpanel::rp.radiogroup(smooth.panel, method,
                      c("aicc", "cv", "manual"), title = "Choice of bandwidth",
                      action = set.bandwidth)
   rpanel::rp.slider(smooth.panel, h.manual, hvec[1], hvec[nvec], plotfun, "df", log = TRUE)
   rpanel::rp.checkbox(smooth.panel, se, plotfun, title = "Standard errors")
   rpanel::rp.radiogroup(smooth.panel, model,
                      c("none", "no effect", "linear"),
                      title = "Reference model", action = plotfun)
   rpanel::rp.checkbox(smooth.panel, test, plotfun, title = "Test")
   rpanel::rp.do(smooth.panel, plotfun)
       
   invisible(smooth.panel)
   }


#' @keywords internal
plott.smooth2 <- function(panel) {
	
   if (panel$method != panel$method.old) {
      if (panel$method != "manual") {
         panel$h <- h.select(panel$x, panel$y, panel$weights, method = panel$method)
         panel$h <- panel$h[1]
       }
      panel$method.old <- panel$method
   }
   
   if (panel$method == "manual")
      panel$h <- panel$h.manual
   if (panel$structure.2d == "common")
      h2 <- rep(panel$h, 2)
   else 
      h2 <- c(panel$h, panel$h * sqrt(wvar(panel$x[,2], panel$weights) /
                                      wvar(panel$x[,1], panel$weights)))
   panel$opt$display <- panel$display
   if (panel$display == "image") {
      opt1 <- panel$opt
      opt1$display <- "slice"
      opt1$add     <- TRUE
      }
   if (panel$display != "rgl")
      panel$surf.ids <- rep(NA, 2)
   panel$opt$theta   <- panel$theta
   panel$opt$phi     <- panel$phi
   panel$opt$se      <- panel$se
   panel$opt$test    <- panel$test
   if (panel$model == "none") panel$opt$test <- FALSE
   if (panel$display == "rgl") {
   	  if (panel$display.old == "rgl")
         panel$opt$add <- TRUE
   	  else
   	     panel$opt$add <- FALSE
      }
   result <- sm.regression.2d(panel$x, panel$y, h2, panel$model, panel$weights, 
                    panel$rawdata, panel$opt)
   if (panel$display == "rgl") {
      if (!is.na(sum(panel$surf.ids)))
         rgl::pop3d(id = panel$surf.ids)
      panel$surf.ids <- result$surf.ids
   }
   if (panel$display == "rgl" & !panel$opt$add) panel$opt$scaling <- result$scaling
   with(panel, {
      if (display == "image")
         sm.regression.2d(x, y, h2, model, weights, rawdata, opt1)
      if (!(display == "rgl")) {
         df.h    <- approx(hvec, dfvec, h, rule = 2)$y
         df.text <- paste("df =", as.character(round(df.h, 1)))
         h.text  <- paste("   h = (", as.character(signif(h2[1], 3)),",",
                                      as.character(signif(h2[2], 3)), ")")
         if ("p" %in% names(result))
            p.text <- paste("   p =", round(result$p, 3))
         else
            p.text <- ""
         title(paste(df.text, h.text, p.text))
      }
      else if ("p" %in% names(result))
         cat("Test (", panel$model, "):  p = ", round(result$p, 3), "\n", sep = "") 
   })
   panel$display.old <- panel$display
   panel
   }
   
replott.smooth2 <- function(panel) {
   rpanel::rp.tkrreplot(panel, smplot)
   panel
   }

rp.smooth2 <- function(x, y, h, model, weights, rawdata, opt) {

   opt <- sm.options(opt)
   replace.na(opt, se,         FALSE)
   replace.na(opt, test,       FALSE)
   if (is.na(opt$display)) display.set <- FALSE
      else                 display.set <- TRUE
   replace.na(opt, display, "persp")
   if (!display.set | (display.set & opt$display == "rgl"))
      opt$panel.plot <- FALSE
   else
      replace.na(opt, panel.plot, TRUE)    
   if (opt$structure.2d == "different")
      stop("structure.2d cannot be set to different when panel = TRUE.")
   opt$verbose <- 0
   nvec  <- 20
   hvec  <- rep(0, nvec)
   hvec[1]    <- h.select(x, y, weights = weights, method = "df", df = 2.1, 
                    structure.2d = opt$structure.2d, nbins = 0)[1]
   hvec[nvec] <- h.select(x, y, weights = weights, method = "df", df = 30,
                    structure.2d = opt$structure.2d, nbins = 0)[1]
   hvec       <- exp(seq(log(hvec[1]), log(hvec[nvec]), length = nvec))
   dfvec      <- seq(2.1, 30, length = nvec)
   for (i in 2:(nvec - 1)) {
   	  if (opt$structure.2d == "common")
         h2 <- rep(hvec[i], 2)
      else 
         h2 <- c(hvec[i], hvec[i] * sqrt(wvar(x[,2], weights) /
                                         wvar(x[,1], weights)))
      dfvec[i] <- sum(diag(sm.weight2(x, x, h2, weights = weights)))
      }
   if (opt$panel.plot & !requireNamespace("tkrplot", quietly = TRUE)) opt$panel.plot <- FALSE
   
   smooth.panel <- rpanel::rp.control("Nonparametric regression - 2 covariates", 
                      x = x, y = y, h = h[1], structure.2d = opt$structure.2d,
                      model = model, weights = weights, rawdata = rawdata,
                      opt = opt, hvec = hvec, dfvec = dfvec, h.manual = h[1], 
                      display = opt$display, display.old = "none", surf.ids = rep(NA, 2),
                      theta = opt$theta, phi = opt$phi,
                      method = "manual", method.old = "manual", se = opt$se, test = opt$test)
   if (opt$panel.plot) {
      rpanel::rp.tkrplot(smooth.panel, smplot, plott.smooth2, pos = "right",
                 hscale = opt$hscale, vscale = opt$vscale)
      plotfun <- replott.smooth2
      }
   else
      plotfun <- plott.smooth2
   rpanel::rp.radiogroup(smooth.panel, method,
                      c("aicc", "cv", "manual"), title = "Choice of bandwidth",
                      action = plotfun)
   rpanel::rp.slider(smooth.panel, h.manual, hvec[1], hvec[nvec], 
                      plotfun, "df", log = TRUE)
   rpanel::rp.checkbox(smooth.panel, se, plotfun, title = "Standard errors")
   rpanel::rp.radiogroup(smooth.panel, model,
                      c("none", "no effect", "linear"),
                      title = "Reference model", action = plotfun)
   rpanel::rp.checkbox(smooth.panel, test, plotfun, title = "Test")
   if (!display.set) {
      display.options <- c("persp", "image")
      if (requireNamespace("rgl", quietly = TRUE)) 
         display.options <- c(display.options, "rgl")
      rpanel::rp.radiogroup(smooth.panel, display, display.options,
                      title = "Display", action = plotfun)
      }
   if (opt$display == "persp") {
      rpanel::rp.slider(smooth.panel, theta, -180, 180, plotfun, "persp angle 1")
      rpanel::rp.slider(smooth.panel, phi,      0,  90, plotfun, "persp angle 2")
      }
       
   invisible(smooth.panel)
   }

#     Density estimation

#' @keywords internal
plott.density1 <- function(panel) {
   if (panel$method == "manual") panel$h <- panel$h.manual
   panel$opt$se   <- panel$se
   if (panel$se & panel$model == "normal") 
      panel$opt$band <- TRUE
   else
      panel$opt$band <- FALSE
   with(panel, {
      result <- sm.density.1d(x, h, model, weights, rawdata, opt)
      h.text  <- paste("   h =", as.character(signif(h, 3)))
      title(h.text)
      })
   panel
   }
   
set.bandwidth.d <- function(panel) {
   if (panel$method != "manual") {
      panel$h <- h.select(panel$x, NA, panel$weights, method = panel$method)
      panel$h <- panel$h[1]
      }
   if (is.matrix(panel$x)) ndim <- ncol(panel$x)
      else                 ndim <- 1
   if (panel$opt$panel.plot) {
      if      (ndim == 1) rpanel::rp.do(panel, replott.density1)
      else if (ndim == 2) rpanel::rp.do(panel, replott.density2)
      else                rpanel::rp.do(panel, replott.density3)
   }
   else {
      if      (ndim == 1) rpanel::rp.do(panel, plott.density1)
      else if (ndim == 2) rpanel::rp.do(panel, plott.density2)
      else                rpanel::rp.do(panel, plott.density3)
      }
   panel
   }

replott.density1 <- function(panel) {
   rpanel::rp.tkrreplot(panel, plot)
   panel
   }

rp.density1 <- function(x, h, model, weights, rawdata, opt) {

   opt <- sm.options(opt)
   replace.na(opt, display,    "lines")
   replace.na(opt, se,         FALSE)
   replace.na(opt, panel.plot, TRUE)
   opt$verbose <- 0
   if (opt$panel.plot & !requireNamespace("tkrplot", quietly = TRUE)) opt$panel <- FALSE

   smooth.panel <- rpanel::rp.control("Density estimation - 1 variable", 
                      x = x, h = h,
                      model = model, weights = weights, rawdata = rawdata,
                      opt = opt, h.manual = h,
                      method = "manual", se = opt$se)
   if (opt$panel.plot) {
      rpanel::rp.tkrplot(smooth.panel, plot, plott.density1, pos = "right",
                 hscale = opt$hscale, vscale = opt$vscale)
      plotfun <- replott.density1
      }
   else
      plotfun <- plott.density1
   rpanel::rp.radiogroup(smooth.panel, method,
                      c("normal", "sj", "cv", "manual"), title = "Choice of bandwidth",
                      action = set.bandwidth.d)
   rpanel::rp.slider(smooth.panel, h.manual, h / 10, h * 10, plotfun, "h", log = TRUE)
   rpanel::rp.checkbox(smooth.panel, se, plotfun, title = "Standard errors")
   rpanel::rp.radiogroup(smooth.panel, model, c("none", "normal"),
                      title = "Reference model", action = plotfun)
   rpanel::rp.do(smooth.panel, plotfun)
       
   invisible(smooth.panel)
   }

#' @keywords internal
plott.density2 <- function(panel) {

   if (panel$method != panel$method.old) {
      if (panel$method != "manual") {
         panel$h <- h.select(panel$x, NA, panel$weights, method = panel$method)
         panel$h <- panel$h[1]
       }
      panel$method.old <- panel$method
   }
   
   if (panel$method == "manual")
      panel$h <- panel$h.manual
   if (panel$structure.2d == "common")
      h2 <- rep(panel$h, 2)
   else 
      h2 <- c(panel$h, panel$h * sqrt(wvar(panel$x[,2], panel$weights) /
                                      wvar(panel$x[,1], panel$weights)))
   panel$opt$display <- panel$display
   if (panel$display == "image") {
      opt1 <- panel$opt
      opt1$display <- "slice"
      opt1$add     <- TRUE
      }
   if (panel$display != "rgl")
      panel$surf.ids <- rep(NA, 2)
   panel$opt$theta   <- panel$theta
   panel$opt$phi     <- panel$phi
   panel$opt$se      <- panel$se
   panel$opt$test    <- panel$test
   if (panel$model == "none") panel$opt$test <- FALSE
   if (panel$display == "rgl") {
   	  if (panel$display.old == "rgl")
         panel$opt$add <- TRUE
   	  else
   	     panel$opt$add <- FALSE
      }
   result <- sm.density.2d(panel$x, h2, panel$weights, panel$rawdata, panel$opt)
   if (panel$display == "rgl") {
      if (!is.na(sum(panel$surf.ids)))
         rgl::pop3d(id = panel$surf.ids)
      panel$surf.ids <- result$surf.ids
   }

   if (panel$display == "rgl" & !panel$opt$add) panel$opt$scaling <- result$scaling
   with(panel, {
      if (display == "image") sm.density.2d(x, h2, weights, rawdata, opt1)
      if (!(display == "rgl")) {
         h.text  <- paste("   h = (", as.character(signif(h2[1], 3)),",",
                                      as.character(signif(h2[2], 3)), ")")
         title(h.text)
         }
      })
   panel$display.old <- panel$display
   panel
   }
   
replott.density2 <- function(panel) {
   rpanel::rp.tkrreplot(panel, smplot)
   panel
   }

rp.density2 <- function(x, h, model, weights, rawdata, opt) {

   opt <- sm.options(opt)
   replace.na(opt, se,         FALSE)
   replace.na(opt, test,       FALSE)
   if (is.na(opt$display)) display.set <- FALSE
      else                 display.set <- TRUE
   replace.na(opt, display, "persp")
   if (!display.set | (display.set & opt$display == "rgl"))
      opt$panel.plot <- FALSE
   else
      replace.na(opt, panel.plot, TRUE)    
   opt$verbose <- 0
   if (opt$panel.plot & !requireNamespace("tkrplot", quietly = TRUE)) opt$panel <- FALSE
   
   smooth.panel <- rpanel::rp.control("Density estimation - 2 variables", 
                      x = x, h = h[1], structure.2d = opt$structure.2d,
                      model = model, weights = weights, rawdata = rawdata,
                      opt = opt, h.manual = h[1], 
                      display = opt$display, display.old = "none",
                      theta = opt$theta, phi = opt$phi, surf.ids = rep(NA, 2),
                      method = "manual", method.old = "manual", se = opt$se, test = opt$test)

   if (opt$panel.plot) {
      rpanel::rp.tkrplot(smooth.panel, smplot, plott.density2, pos = "right",
                 hscale = opt$hscale, vscale = opt$vscale)
      plotfun <- replott.density2
      }
   else
      plotfun <- plott.density2
   rpanel::rp.radiogroup(smooth.panel, method,
                      c("normal", "cv", "manual"), title = "Choice of bandwidth",
                      action = plotfun)
   rpanel::rp.slider(smooth.panel, h.manual, h[1] / 10, h[1] * 10, 
                      plotfun, "h", log = TRUE)
   # rpanel::rp.checkbox(smooth.panel, se, plotfun, title = "Standard errors")
   # rpanel::rp.radiogroup(smooth.panel, model, c("none", "normal"),
   #                    title = "Reference model", action = plotfun)
   if (!display.set) {
      display.options <- c("persp", "image")
      if (requireNamespace("rgl", quietly = TRUE)) 
         display.options <- c(display.options, "rgl")
      rpanel::rp.radiogroup(smooth.panel, display, display.options,
                      title = "Display", action = plotfun)
      }
   if (opt$display == "persp") {
      rpanel::rp.slider(smooth.panel, theta, -180, 180, plotfun, "persp angle 1")
      rpanel::rp.slider(smooth.panel, phi,      0,  90, plotfun, "persp angle 2")
      }
       
   invisible(smooth.panel)
   }

#' @keywords internal
plott.density3 <- function(panel) {
   if (panel$method != panel$method.old) {
      if (panel$method != "manual") {
         panel$h <- h.select(panel$x, panel$y, panel$weights, method = panel$method)
         panel$h <- panel$h[1]
       }
      panel$method.old <- panel$method
   }

   if (panel$method == "manual") panel$h <- panel$h.manual
   if (panel$structure.2d == "common")
      h3 <- rep(panel$h, 3)
   else 
      h3 <- c(panel$h, panel$h * sqrt(wvar(panel$x[,2], panel$weights) /
                                      wvar(panel$x[,1], panel$weights)),
                       panel$h * sqrt(wvar(panel$x[,3], panel$weights) /
                                      wvar(panel$x[,1], panel$weights)))
   panel$opt$display <- panel$display
   if (panel$display.old == "rgl")
      panel$opt$add <- TRUE
   else
   	  panel$opt$add <- FALSE
   result <- sm.density.3d(panel$x, h3, panel$weights, panel$rawdata, panel$opt)
   if (panel$display == "rgl" & !panel$opt$add) panel$opt$scaling <- result$scaling
   if (!all(is.na(panel$surf.ids)))
      rgl::pop3d(id = panel$surf.ids)
   panel$surf.ids <- result$surf.ids

   panel$display.old <- panel$display
   panel
   }
   
replott.density3 <- function(panel) {
   rpanel::rp.tkrreplot(panel, smplot)
   panel
   }

rp.density3 <- function(x, h, model, weights, rawdata, opt) {

   opt <- sm.options(opt)
   replace.na(opt, se,         FALSE)
   replace.na(opt, test,       FALSE)
   if (is.na(opt$display)) display.set <- FALSE
      else                 display.set <- TRUE
   replace.na(opt, display, "rgl")
   opt$panel.plot <- FALSE
   opt$verbose <- 0
   
   smooth.panel <- rpanel::rp.control("Density estimation - 3 variables", 
                      x = x, h = h[1], structure.2d = opt$structure.2d,
                      model = model, weights = weights, rawdata = rawdata,
                      opt = opt, h.manual = h[1], 
                      display = opt$display, display.old = "none",
                      theta = opt$theta, phi = opt$phi, surf.ids = NA,
                      method = "manual", method.old = "manual", se = opt$se, test = opt$test)
   if (opt$panel.plot) {
      rpanel::rp.tkrplot(smooth.panel, smplot, plott.density2, pos = "right",
                 hscale = opt$hscale, vscale = opt$vscale)
      plotfun <- replott.density3
      }
   else
      plotfun <- plott.density3
   rpanel::rp.radiogroup(smooth.panel, method,
                      c("normal", "manual"), title = "Choice of bandwidth",
                      action = plotfun)
   rpanel::rp.slider(smooth.panel, h.manual, h[1] / 3, h[1] * 3, 
                      plotfun, "h", log = TRUE)
       
   invisible(smooth.panel)
   }

