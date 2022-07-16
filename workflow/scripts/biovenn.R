#' This is a significant modification of BioVenn::draw.venn customised for the
#' GRAVI workflow. This may be integrated within extraChIPs at some point.
#' The input `data` shuold be a list. If length(data) > 3, only the 1st three
#' elements will be used
#' @import plotrix
biovenn <- function (
    data,
    title = c(),
    fontfamily = "serif",
    fontcol = "black",
    title.font = 2,
    label.font = 2,
    nbr.font = 1,
    title.cex = 2,
    labels.cex = 1.5,
    nbr.cex = 1.5,
    cols = c("red", "green", "blue"),
    bg_c = "white",
    shift_lab = 1.2
    ) {

  stopifnot(is(data, "list"))
  stopifnot(length(data) > 1)
  list_x <- unique(data[[1]])
  list_y <- unique(data[[2]])
  xtitle <- names(data)[[1]]
  ytitle <- names(data)[[2]]
  list_z <- NULL
  if (length(x) == 3) {
    list_z <- unique(data[[3]])
    ztitle <- names(data)[[3]]
  }
  x <- length(list_x)
  y <- length(list_y)
  z <- length(list_z)
  list_xy <- intersect(list_x, list_y)
  xy <- length(list_xy)
  list_xz <- intersect(list_x, list_z)
  xz <- length(list_xz)
  list_yz <- intersect(list_y, list_z)
  yz <- length(list_yz)
  list_xyz <- intersect(list_xy, list_z)
  xyz <- length(list_xyz)
  list_xy_only <- setdiff(list_xy, list_xyz)
  xy_only <- length(list_xy_only)
  list_xz_only <- setdiff(list_xz, list_xyz)
  xz_only <- length(list_xz_only)
  list_yz_only <- setdiff(list_yz, list_xyz)
  yz_only <- length(list_yz_only)
  list_x_only <- setdiff(list_x, c(list_xy, list_xz))
  x_only <- length(list_x_only)
  list_y_only <- setdiff(list_y, c(list_xy, list_yz))
  y_only <- length(list_y_only)
  list_z_only <- setdiff(list_z, c(list_xz, list_yz))
  z_only <- length(list_z_only)

  cols <- rep(cols, 3)
  x_c <- cols[[1]]
  y_c <- cols[[2]]
  z_c <- cols[[3]]

  sq <- function(nr) nr^2
  sqr <- function(nr) sqrt(round(abs(nr)))
  arccos <- function(nr) {
    acos(min(max(-1, round(nr, 5)), 1))
  }
  width_p = 1000
  height_p = 1000
  amp = 1e+05/(x + y + z - xy - xz - yz + xyz)
  x_text = x
  x = x * amp
  y_text = y
  y = y * amp
  z_text = z
  z = z * amp
  xy_text = xy
  xy = xy * amp
  xz_text = xz
  xz = xz * amp
  yz_text = yz
  yz = yz * amp
  xyz_text = xyz
  xyz = xyz * amp
  total = x + y + z - xy - xz - yz + xyz
  total_text = x_text + y_text + z_text - xy_text - xz_text -
    yz_text + xyz_text
  x_r = sqr(x/pi)
  y_r = sqr(y/pi)
  z_r = sqr(z/pi)
  xy_d = x_r + y_r
  if (x && y) {
    while (
      xy * 0.999 > sq(x_r) * arccos((sq(xy_d) + sq(x_r) -
                                     sq(y_r))/(2 * xy_d * x_r)) + sq(y_r) * arccos((sq(xy_d) +
                                                                                    sq(y_r) - sq(x_r))/(2 * xy_d * y_r)) - 0.5 * sqr(round((xy_d +
                                                                                                                                            x_r + y_r) * (xy_d + x_r - y_r) * (xy_d - x_r + y_r) *
                                                                                                                                           (-xy_d + x_r + y_r), 5))) {
      xy_d = xy_d - min(x_r, y_r)/1000
    }
  }
  xz_d = x_r + z_r
  if (x && z) {
    while (xz * 0.999 > sq(x_r) * arccos((sq(xz_d) + sq(x_r) -
                                          sq(z_r))/(2 * xz_d * x_r)) + sq(z_r) * arccos((sq(xz_d) +
                                                                                         sq(z_r) - sq(x_r))/(2 * xz_d * z_r)) - 0.5 * sqr(round((xz_d +
                                                                                                                                                 x_r + z_r) * (xz_d + x_r - z_r) * (xz_d - x_r + z_r) *
                                                                                                                                                (-xz_d + x_r + z_r), 5))) {
      xz_d = xz_d - min(x_r, z_r)/1000
    }
  }
  yz_d = y_r + z_r
  if (y && z) {
    while (yz * 0.999 > sq(y_r) * arccos((sq(yz_d) + sq(y_r) -
                                          sq(z_r))/(2 * yz_d * y_r)) + sq(z_r) * arccos((sq(yz_d) +
                                                                                         sq(z_r) - sq(y_r))/(2 * yz_d * z_r)) - 0.5 * sqr(round((yz_d +
                                                                                                                                                 y_r + z_r) * (yz_d + y_r - z_r) * (yz_d - y_r + z_r) *
                                                                                                                                                (-yz_d + y_r + z_r), 5))) {
      yz_d = yz_d - min(y_r, z_r)/1000
    }
  }
  if (xy_d > xz_d + yz_d) {
    xy_d = xz_d + yz_d
  }
  if (xz_d > xy_d + yz_d) {
    xz_d = xy_d + yz_d
  }
  if (yz_d > xy_d + xz_d) {
    yz_d = xy_d + xz_d
  }
  x_a = arccos((sq(xy_d) + sq(xz_d) - sq(yz_d))/(2 * xy_d * xz_d))
  y_a = arccos((sq(xy_d) + sq(yz_d) - sq(xz_d))/(2 * xy_d * yz_d))
  z_a = arccos((sq(xz_d) + sq(yz_d) - sq(xy_d))/(2 * xz_d * yz_d))
  x_yz = xz_d * sin(z_a)
  y_yz = xy_d * cos(y_a)
  width_h = max(y_r + y_yz, x_r, z_r - yz_d + y_yz) + max(x_r, y_r - y_yz, z_r + yz_d - y_yz)
  ppu_h = width_p/width_h
  width_v = max(x_r + x_yz, y_r, z_r) + max(y_r, z_r, x_r - x_yz)
  ppu_v = height_p/width_v
  ppu = min(ppu_h, ppu_v)
  x_h = max(x_r, y_r + y_yz, z_r - yz_d + y_yz)
  x_v = max(x_r, y_r - x_yz, z_r - x_yz)
  y_h = max(x_r - y_yz, y_r, z_r - yz_d)
  y_v = max(x_r + x_yz, y_r, z_r)
  z_h = max(x_r + yz_d - y_yz, y_r + yz_d, z_r)
  z_v = max(x_r + x_yz, y_r, z_r)

  xy_i_h_part1 = (x_h + y_h)/2 + ((y_h - x_h) * (sq(x_r) - sq(y_r)))/(2 * sq(xy_d))
  xy_i_v_part1 = (x_v + y_v)/2 + ((y_v - x_v) * (sq(x_r) - sq(y_r)))/(2 * sq(xy_d))
  xy_i_h_part2 = 2 * ((x_v - y_v)/sq(xy_d)) * sqr((xy_d + x_r + y_r) * (xy_d + x_r - y_r) * (xy_d - x_r + y_r) * (-xy_d + x_r + y_r))/4
  xy_i_v_part2 = 2 * ((x_h - y_h)/sq(xy_d)) * sqr((xy_d + x_r + y_r) * (xy_d + x_r - y_r) * (xy_d - x_r + y_r) * (-xy_d + x_r + y_r))/4

  xy_i1_h = xy_i_h_part1 - xy_i_h_part2
  xy_i1_v = xy_i_v_part1 + xy_i_v_part2
  xy_i2_h = xy_i_h_part1 + xy_i_h_part2
  xy_i2_v = xy_i_v_part1 - xy_i_v_part2

  xz_i_h_part1 = (x_h + z_h)/2 + ((z_h - x_h) * (sq(x_r) - sq(z_r)))/(2 * sq(xz_d))
  xz_i_v_part1 = (x_v + z_v)/2 + ((z_v - x_v) * (sq(x_r) - sq(z_r)))/(2 * sq(xz_d))
  xz_i_h_part2 = 2 * ((x_v - z_v)/sq(xz_d)) * sqr((xz_d + x_r + z_r) * (xz_d + x_r - z_r) * (xz_d - x_r + z_r) * (-xz_d + x_r + z_r))/4
  xz_i_v_part2 = 2 * ((x_h - z_h)/sq(xz_d)) * sqr((xz_d + x_r + z_r) * (xz_d + x_r - z_r) * (xz_d - x_r + z_r) * (-xz_d + x_r + z_r))/4

  xz_i1_h = xz_i_h_part1 + xz_i_h_part2
  xz_i1_v = xz_i_v_part1 - xz_i_v_part2
  xz_i2_h = xz_i_h_part1 - xz_i_h_part2
  xz_i2_v = xz_i_v_part1 + xz_i_v_part2

  yz_i_h_part1 = (y_h + z_h)/2 + ((z_h - y_h) * (sq(y_r) - sq(z_r)))/(2 * sq(yz_d))
  yz_i_v_part1 = (y_v + z_v)/2 + ((z_v - y_v) * (sq(y_r) - sq(z_r)))/(2 * sq(yz_d))
  yz_i_h_part2 = 2 * ((y_v - z_v)/sq(yz_d)) * sqr((yz_d + y_r + z_r) * (yz_d + y_r - z_r) * (yz_d - y_r + z_r) * (-yz_d + y_r + z_r))/4
  yz_i_v_part2 = 2 * ((y_h - z_h)/sq(yz_d)) * sqr((yz_d + y_r + z_r) * (yz_d + y_r - z_r) * (yz_d - y_r + z_r) * (-yz_d + y_r + z_r))/4

  yz_i1_h = yz_i_h_part1 - yz_i_h_part2
  yz_i1_v = yz_i_v_part1 + yz_i_v_part2
  yz_i2_h = yz_i_h_part1 + yz_i_h_part2
  yz_i2_v = yz_i_v_part1 - yz_i_v_part2
  if (x && y && z) {
    xyz_f_h = (xy_i1_h + xz_i1_h + yz_i1_h)/3
    xyz_f_v = (xy_i1_v + xz_i1_v + yz_i1_v)/3
  }
  if (x && y && z && xy && xz) {
    xyz_yz_i1 = sqr(sq(xyz_f_h - yz_i1_h) + sq(xyz_f_v - yz_i1_v))
    x_ratio_h = (xyz_f_h - yz_i1_h)/xyz_yz_i1
    x_ratio_v = (xyz_f_v - yz_i1_v)/xyz_yz_i1
    x_out_h = x_h - x_r * x_ratio_h
    x_out_v = x_v - x_r * x_ratio_v
    x_f_h = (x_out_h + yz_i1_h)/2
    x_f_v = (x_out_v + yz_i1_v)/2
  }
  else if (x && y && !z || x && y && z && !xz) {
    xy_f_h = (xy_i1_h + xy_i2_h)/2
    xy_f_v = (xy_i1_v + xy_i2_v)/2
    x_in_h = y_h + cos(y_a) * y_r
    x_in_v = y_v - sin(y_a) * y_r
    x_out_h = x_h + cos(y_a) * x_r
    x_out_v = x_v - sin(y_a) * x_r
    x_f_h = (x_out_h + x_in_h)/2
    x_f_v = (x_out_v + x_in_v)/2
  }
  else if (x && !y && z || x && y && z && !xy) {
    xz_f_h = (xz_i1_h + xz_i2_h)/2
    xz_f_v = (xz_i1_v + xz_i2_v)/2
    x_in_h = z_h - cos(z_a) * z_r
    x_in_v = z_v - sin(z_a) * z_r
    x_out_h = x_h - cos(z_a) * x_r
    x_out_v = x_v - sin(z_a) * x_r
    x_f_h = (x_out_h + x_in_h)/2
    x_f_v = (x_out_v + x_in_v)/2
  }
  if (x && y && z && xy && yz) {
    xyz_xz_i1 = sqr(sq(xyz_f_h - xz_i1_h) + sq(xyz_f_v - xz_i1_v))
    y_ratio_h = (xyz_f_h - xz_i1_h)/xyz_xz_i1
    y_ratio_v = (xyz_f_v - xz_i1_v)/xyz_xz_i1
    y_out_h = y_h - y_r * y_ratio_h
    y_out_v = y_v - y_r * y_ratio_v
    y_f_h = (y_out_h + xz_i1_h)/2
    y_f_v = (y_out_v + xz_i1_v)/2
  }
  else if (x && y && !z || x && y && z && !yz) {
    xy_f_h = (xy_i1_h + xy_i2_h)/2
    xy_f_v = (xy_i1_v + xy_i2_v)/2
    y_in_h = x_h - cos(y_a) * x_r
    y_in_v = x_v + sin(y_a) * x_r
    y_out_h = y_h - cos(y_a) * y_r
    y_out_v = y_v + sin(y_a) * y_r
    y_f_h = (y_out_h + y_in_h)/2
    y_f_v = (y_out_v + y_in_v)/2
  }
  else if (!x && y && z || x && y && z && !xy) {
    yz_f_h = (yz_i1_h + yz_i2_h)/2
    yz_f_v = (yz_i1_v + yz_i2_v)/2
    y_in_h = z_h - z_r
    y_in_v = z_v
    y_out_h = y_h - y_r
    y_out_v = y_v
    y_f_h = (y_out_h + y_in_h)/2
    y_f_v = (y_out_v + y_in_v)/2
  }
  if (x && y && z && xz && yz) {
    xyz_xy_i1 = sqr(sq(xyz_f_h - xy_i1_h) + sq(xyz_f_v - xy_i1_v))
    z_ratio_h = (xyz_f_h - xy_i1_h)/xyz_xy_i1
    z_ratio_v = (xyz_f_v - xy_i1_v)/xyz_xy_i1
    z_out_h = z_h - z_r * z_ratio_h
    z_out_v = z_v - z_r * z_ratio_v
    z_f_h = (z_out_h + xy_i1_h)/2
    z_f_v = (z_out_v + xy_i1_v)/2
  }
  else if (x && !y && z || x && y && z && !yz) {
    xz_f_h = (xz_i1_h + xz_i2_h)/2
    xz_f_v = (xz_i1_v + xz_i2_v)/2
    z_in_h = x_h + cos(z_a) * x_r
    z_in_v = x_v + sin(z_a) * x_r
    z_out_h = z_h + cos(z_a) * z_r
    z_out_v = z_v + sin(z_a) * z_r
    z_f_h = (z_out_h + z_in_h)/2
    z_f_v = (z_out_v + z_in_v)/2
  }
  else if (!x && y && z || x && y && z && !xz) {
    yz_f_h = (yz_i1_h + yz_i2_h)/2
    yz_f_v = (yz_i1_v + yz_i2_v)/2
    z_in_h = y_h + z_r
    z_in_v = y_v
    z_out_h = z_h + y_r
    z_out_v = z_v
    z_f_h = (z_out_h + z_in_h)/2
    z_f_v = (z_out_v + z_in_v)/2
  }
  if (x && y && z) {
    dh = (xyz_f_h - z_h) - (xy_i2_h - z_h)
    dv = (xyz_f_v - z_v) - (xy_i2_v - z_v)
    dr = sqr(sq(dh) + sq(dv))
    D = (xy_i2_h - z_h) * (xyz_f_v - z_v) - (xyz_f_h - z_h) * (xy_i2_v - z_v)
    z_in_h = z_h + (D * dv - dh * sqr(sq(z_r) * sq(dr) - sq(D)))/sq(dr)
    z_in_v = z_v + (-D * dh - abs(dv) * sqr(sq(z_r) * sq(dr) - sq(D)))/sq(dr)
    xy_f_h = (z_in_h + xy_i2_h)/2
    xy_f_v = (z_in_v + xy_i2_v)/2
  }
  if (x && y && z) {
    dh = (xyz_f_h - y_h) - (xz_i2_h - y_h)
    dv = (xyz_f_v - y_v) - (xz_i2_v - y_v)
    dr = sqr(sq(dh) + sq(dv))
    D = (xz_i2_h - y_h) * (xyz_f_v - y_v) - (xyz_f_h - y_h) * (xz_i2_v - y_v)
    y_in_h = y_h + (D * dv - dh * sqr(sq(y_r) * sq(dr) - sq(D)))/sq(dr)
    y_in_v = y_v + (-D * dh - abs(dv) * sqr(sq(y_r) * sq(dr) - sq(D)))/sq(dr)
    xz_f_h = (y_in_h + xz_i2_h)/2
    xz_f_v = (y_in_v + xz_i2_v)/2
  }
  if (x && y && z) {
    dh = (xyz_f_h - x_h) - (yz_i2_h - x_h)
    dv = (xyz_f_v - x_v) - (yz_i2_v - x_v)
    dr = sqr(sq(dh) + sq(dv))
    D = (yz_i2_h - x_h) * (xyz_f_v - x_v) - (xyz_f_h - x_h) * (yz_i2_v - x_v)
    x_in_h = x_h + (D * dv - dh * sqr(sq(x_r) * sq(dr) - sq(D)))/sq(dr)
    x_in_v = x_v + (-D * dh + abs(dv) * sqr(sq(x_r) * sq(dr) - sq(D)))/sq(dr)
    yz_f_h = (x_in_h + yz_i2_h)/2
    yz_f_v = (x_in_v + yz_i2_v)/2
  }
  if (xy_d == xz_d + yz_d || xz_d == xy_d + yz_d || yz_d == xy_d + xz_d) {
    if (x && !x_only && y && !y_only) {
      xz_f_v = yz_f_v = xyz_f_v = x_v
      xz_f_h = (max(y_h + y_r, x_h - x_r) + (x_h + x_r))/2
      yz_f_h = ((y_h - y_r) + min(y_h + y_r, x_h - x_r))/2
    }
    else if (x && !x_only && z && !z_only) {
      xy_f_v = yz_f_v = xyz_f_v = x_v
      xy_f_h = ((x_h - x_r) + min(x_h + x_r, z_h - z_r))/2
      yz_f_h = (max(x_h + x_r, z_h - z_r) + (z_h + z_r))/2
    }
    else if (y && !y_only && z && !z_only) {
      yz_f_v = xz_f_v = xyz_f_v = x_v
      yz_f_h = (max(x_h + x_r, y_h - y_r) + (y_h + y_r))/2
      xz_f_h = (max(y_h + y_r, z_h - z_r) + (z_h + z_r))/2
    }
    else if (x && !x_only) {
      yz_f_v = xyz_f_v = x_v
      if (!xz_only) {
        z_f_h = (max(x_h + x_r, y_h + y_r) + (z_h + z_r))/2
        z_f_v = x_v
        xy_f_h = (max(y_h - y_r, x_h - x_r) + (z_h - z_r))/2
        xy_f_v = x_v
        yz_f_h = (max(x_h + x_r, z_h - z_r) + (y_h + y_r))/2
      }
      else if (!xy_only) {
        y_f_h = ((y_h - y_r) + min(z_h - z_r, x_h - x_r))/2
        y_f_v = x_v
        xz_f_h = (max(y_h + y_r, x_h - x_r) + (x_h + x_r))/2
        xz_f_v = x_v
        yz_f_h = ((z_h - z_r) + min(y_h + y_r, x_h - x_r))/2
      }
    }
    else if (y && !y_only) {
      xz_f_v = xyz_f_v = x_v
      if (!yz_only) {
        z_f_h = (max(y_h + y_r, x_h + x_r) + (z_h + z_r))/2
        z_f_v = x_v
        xy_f_h = ((y_h - y_r) + min(y_h + y_r, z_h - z_r))/2
        xy_f_v = x_v
        xz_f_h = (max(y_h + y_r, z_h - z_r) + (x_h + x_r))/2
      }
      else if (!xy_only) {
        x_f_h = (max(y_h + y_r, z_h + z_r) + (x_h + x_r))/2
        x_f_v = x_v
        yz_f_h = ((y_h - y_r) + max(y_h + y_r, x_h - x_r))/2
        yz_f_v = x_v
        xz_f_h = (max(y_h + y_r, x_h - x_r) + (z_h + z_r))/2
      }
    }
    else if (z && !z_only) {
      xy_f_v = xyz_f_v = x_v
      if (!yz_only) {
        y_f_h = ((y_h - y_r) + min(x_h - x_r, z_h - z_r))/2
        y_f_v = x_v
        xz_f_h = (max(y_h + y_r, z_h - z_r) + (z_h + z_r))/2
        xz_f_v = x_v
        xy_f_h = ((x_h - x_r) + min(y_h + y_r, z_h - z_r))/2
      }
      else if (!xz_only) {
        x_f_h = ((x_h - x_r) + min(y_h - y_r, z_h - z_r))/2
        x_f_v = x_v
        yz_f_h = (max(x_h + x_r, z_h - z_r) + (z_h + z_r))/2
        yz_f_v = x_v
        xy_f_h = ((y_h - y_r) + min(x_h + x_r, z_h - z_r))/2
      }
    }
    xyz_f_h = 0.5 * (
      max(x_h - x_r, y_h - y_r, z_h - z_r) + min(x_h + x_r, y_h + y_r, z_h + z_r)
    )
  }
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mar = c(1, 1, 4, 1))
  graphics::par(pty = "s", bg = bg_c)
  graphics::plot(
    0, type = "n",
    axes = FALSE, xlim = c(0, width_p), ylim = c(height_p, 0),
    xlab = "", ylab = "", xaxt = "none", yaxt = "none"
  )
  graphics::par(family = fontfamily)
  graphics::title(
    main = title, line = 1, font.main = title.font, cex.main = title.cex, col.main = fontcol
  )
  plotrix::draw.circle(
    x = ppu * x_h,
    y = ppu * x_v,
    radius = ppu * x_r,
    border = "black", lty = 1,
    col = grDevices::rgb(
      grDevices::col2rgb(x_c)[, 1][1],
      grDevices::col2rgb(x_c)[, 1][2],
      grDevices::col2rgb(x_c)[, 1][3],
      maxColorValue = 255, alpha = 128)
  )
  plotrix::draw.circle(
    ppu * y_h, ppu * y_v, ppu * y_r,
    border = "black",lty = 1,
    col = grDevices::rgb(
      grDevices::col2rgb(y_c)[, 1][1],
      grDevices::col2rgb(y_c)[, 1][2],
      grDevices::col2rgb(y_c)[, 1][3],
      maxColorValue = 255, alpha = 128)
  )
  plotrix::draw.circle(
    ppu * z_h, ppu * z_v, ppu * z_r,
    border = "black",lty = 1,
    col = grDevices::rgb(
      grDevices::col2rgb(z_c)[, 1][1],
      grDevices::col2rgb(z_c)[, 1][2],
      grDevices::col2rgb(z_c)[, 1][3],
      maxColorValue = 255, alpha = 128)
  )
  if (x_only) {
    graphics::text(
      ppu * x_f_h, ppu * x_f_v, x_only,
      adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (y_only) {
    graphics::text(
      ppu * y_f_h, ppu * y_f_v, y_only,
      adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (z_only) {
    graphics::text(
      ppu * z_f_h, ppu * z_f_v, z_only,
      adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (xy_only) {
    graphics::text(
      ppu * xy_f_h, ppu * xy_f_v, xy_only,
      adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (xz_only) {
    graphics::text(
      ppu * xz_f_h, ppu * xz_f_v, xz_only,
      adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (yz_only) {
    graphics::text(
      ppu * yz_f_h, ppu * yz_f_v, yz_only,
      adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (xyz) {
    graphics::text(
      ppu * xyz_f_h, ppu * xyz_f_v,
      xyz_text, adj = c(0.5, 0.5), col = fontcol, family = fontfamily,
      font = nbr.font, cex = nbr.cex
    )
  }
  if (x) {
    graphics::text(
      x = ppu * x_h,
      y = ppu * x_v - ppu * x_r * shift_lab / 1.15,
      xtitle,
      adj = c(0.5, 0),
      col = fontcol, family = fontfamily, font = label.font, cex = labels.cex
    )
  }
  if (y) {
    graphics::text(
      x = ppu * y_h - shift_lab * 2 * ppu * y_r / 3,
      y = ppu * y_v + shift_lab * 2 * ppu * y_r / 3,
      ytitle,
      adj = c(0, 0),
      col = fontcol, family = fontfamily, font = label.font, cex = labels.cex
    )
  }
  if (z) {
    graphics::text(
      x = ppu * z_h + shift_lab * 2 * ppu * z_r / 3,
      y = ppu * z_v + shift_lab * 2 * ppu * z_r / 3,
      ztitle,
      adj = c(1, 0),
      col = fontcol, family = fontfamily, font = label.font, cex = labels.cex
    )
  }
  invisible(NULL)
}
