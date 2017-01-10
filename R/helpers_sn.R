getheight.sn = function(scale, shape) {
  # See /inst/development/sharpness

  (0.4/scale) * (-4.2702E-5 * shape^4  +  2.1909E-3 * shape^3  -  4.0631E-2 * shape^2  +  3.2862E-1 * shape^1  +  9.332844E-1)

  }

getwidth95th.sn = function(scale) {

  3.28970725390297 * scale

  }

getmean.sn = function(location, scale, shape) {

  d = shape / ( 1 + shape^2)^0.5

  location + scale * d * (2/pi)^0.5

  }
