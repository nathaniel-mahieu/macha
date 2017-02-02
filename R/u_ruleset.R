atm = list(
  e = list(m = 0.000548579909),
  p = list(m = 1.00727647),
  n = list(m = 1.008671087),

  n14 = list(m = 14.0030740048, a = 0.99632),
  n15 = list(m = 15.0001088982, a = 0.003268),

  c12 = list(m = 12.0000, a = 0.9893),
  c13 = list(m = 13.0033548378, a = 0.0107),
  c14 = list(m = 14.003241989, a = 0),

  h1 = list(m = 1.00782503207, a = 0.999885),
  h2 = list(m = 2.0141017778, a = 0.000115),

  o16 = list(m = 15.994915, a = 0.99757),
  o17 = list(m = 16.999132, a = 0.00038),
  o18 = list(m = 17.999160, a = 0.00205),

  li6 = list(m = 6.015122, a = 0.0759),
  li7 = list(m = 7.016004, a = 0.9241),

  na23 = list(m = 22.989770, a = 1),

  k39 = list(m = 38.963707, a = 0.932581),
  k41 = list(m = 40.961826, a = 0.067302),

  cl35 = list(m = 34.968853, a = 0.7578),
  cl37 = list(m = 36.965903, a = 0.2422),

  br79 = list(m = 78.918338, a = 0.5069),
  br81 = list(m = 80.916291, a = 0.4931),

  mg24 = list(m = 23.985042, a = 0.7899),
  mg25 = list(m = 24.985837, a = 0.1000),
  mg26 = list(m = 25.982593, a = 0.1101),

  s32 = list(m = 31.972071, a = 0.9493),
  s33 = list(m = 32.971458, a = 0.0076),
  s34 = list(m = 33.967867, a = 0.0429),

  f19 = list(m = 18.998403, a = 1),

  p31 = list(m = 30.973762, a = 1),

  si28 = list(m = 27.976927, a = 0.922297),
  si29 = list(m = 28.976495, a = 0.046832),
  si30 = list(m = 29.973770, a = 0.030872)
)

makeruletable = function(l, tag) {
  lapply(seq_along(l), function(i) {
    l = l[[i]]
    l$gain$dir = rep("gain", nrow(l$gain))
    l$loss$dir = rep("loss", nrow(l$loss))
    cbind(rbind(l$gain, l$loss), rule = paste(tag, i))
  }) %>% do.call(what = rbind)
}



isotopes = list(
  list(gain = data.frame(name = "n15", m = atm$n15$m, z = 0), loss = data.frame(name = "n14", m = atm$n14$m, z = 0)),
  list(gain = data.frame(name = "o18", m = atm$o18$m, z = 0), loss = data.frame(name = "o16", m = atm$o16$m, z = 0)),
  list(gain = data.frame(name = "s33", m = atm$s33$m, z = 0), loss = data.frame(name = "s32", m = atm$s32$m, z = 0)),
  list(gain = data.frame(name = "s34", m = atm$s34$m, z = 0), loss = data.frame(name = "s32", m = atm$s32$m, z = 0)),
  list(gain = data.frame(name = "cl37", m = atm$cl37$m, z = 0), loss = data.frame(name = "cl35", m = atm$cl35$m, z = 0)),
  list(gain = data.frame(name = "c13", m = atm$c13$m, z = 0), loss = data.frame(name = "c12", m = atm$c12$m, z = 0)),
  list(gain = data.frame(name = "br81", m = atm$br81$m, z = 0), loss = data.frame(name = "br79", m = atm$br79$m, z = 0)),
  list(gain = data.frame(name = "si29", m = atm$si29$m, z = 0), loss = data.frame(name = "si28", m = atm$si28$m, z = 0)),
  list(gain = data.frame(name = "si30", m = atm$si30$m, z = 0), loss = data.frame(name = "si28", m = atm$si28$m, z = 0)),
  list(gain = data.frame(name = "k41", m = atm$k41$m, z = 0), loss = data.frame(name = "k39", m = atm$k39$m, z = 0))
  )


makeruletable(isotopes, "iso")




# Gain of a charge carrier (z=1 -> z=2)
z1 = list(
  list(gain = data.frame(name = "H+", m = atm$h1$m, z = 1), loss = data.frame()),
  list(gain = data.frame(name = "Na+", m = atm$na23$m, z = 1), loss = data.frame()),
  list(gain = data.frame(name = "K+", m = atm$k39$m, z = 1), loss = data.frame()),
  list(gain = data.frame(name = "Cl-", m = atm$cl35$m, z = -1), loss = data.frame()),
  list(gain = data.frame(name = "Br-", m = atm$br79$m, z = -1), loss = data.frame())
)

makeruletable(z1, "z1")



# Change of a charge carrier (Na -> H) (Charge Neutral)
z2 = list(
  list(gain = data.frame(name = "H+", m = atm$h1$m, z = 1), loss = data.frame(name = "Na+", m = atm$na23$m, z = 1)),
  list(gain = data.frame(name = "H+", m = atm$h1$m, z = 1), loss = data.frame(name = "K+", m = atm$k39$m, z = 1)),
  list(gain = data.frame(name = "K+", m = atm$k39$m, z = 1), loss = data.frame(name = "Na+", m = atm$na23$m, z = 1)),

  list(gain = data.frame(name = "NH4+", m = atm$n14$m + atm$h1$m*4, z = 1), loss = data.frame(name = "Na+", m = atm$na23$m, z = 1)), # NH4 necessary I think
  list(gain = data.frame(name = "NH4+", m = atm$n14$m + atm$h1$m*4, z = 1), loss = data.frame(name = "K+", m = atm$k39$m, z = 1)), # NH4 salts necessary I think

  list(gain = data.frame(name = "Cl-", m = atm$cl35$m, z = -1), loss = data.frame(name = "Br-", m = atm$h1$m, z = -1))
)

makeruletable(z2, "z2")


# Change of a charge carrier (Na -> H) (Charge Reversal) (dz = 2)
z3 = list(
  list(gain = rbind(data.frame(name = "H+", m = atm$h1$m, z = 1), data.frame(name = "H+", m = atm$h1$m, z = 1)), loss = data.frame()),
  list(gain = rbind(data.frame(name = "H+", m = atm$h1$m, z = 1), data.frame(name = "K+", m = atm$k39$m, z = 1)), loss = data.frame()),
  list(gain = rbind(data.frame(name = "H+", m = atm$h1$m, z = 1), data.frame(name = "Na+", m = atm$na23$m, z = 1)), loss = data.frame()),

  list(gain = rbind(data.frame(name = "H+", m = atm$h1$m, z = 1)), loss = data.frame(name = "Cl-", m = atm$cl35$m, z = -1)),
  list(gain = rbind(data.frame(name = "H+", m = atm$h1$m, z = 1)), loss = data.frame(name = "Br-", m = atm$br79$m, z = -1)),

  list(gain = rbind(data.frame(name = "Na+", m = atm$na23$m, z = 1)), loss = data.frame(name = "Cl-", m = atm$cl35$m, z = -1)),
  list(gain = rbind(data.frame(name = "Na+", m = atm$na23$m, z = 1)), loss = data.frame(name = "Br-", m = atm$br79$m, z = -1)),

  list(gain = rbind(data.frame(name = "K+", m = atm$k39$m, z = 1)), loss = data.frame(name = "Cl-", m = atm$cl35$m, z = -1)),
  list(gain = rbind(data.frame(name = "K+", m = atm$k39$m, z = 1)), loss = data.frame(name = "Br-", m = atm$br79$m, z = -1))
)


makeruletable(z3, "z3")

n1 = list(
  list(gain = data.frame(), loss = data.frame(name = "H2O", m = atm$h1$m*2+atm$o16$m, z = 0)),
  list(gain = data.frame(), loss = data.frame(name = "CO2", m = atm$c12$m + atm$o16$m*2, z = 0)),
  list(gain = data.frame(), loss = data.frame(name = "NH3", m = atm$n14$m + atm$h1$m*3, z = 0)),
  list(gain = data.frame(), loss = data.frame(name = "CO", m = atm$c12$m + atm$o16$m, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "HCOOH", m = atm$c12$m*2 + atm$h1$m*4 + atm$o16$m*2, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "CH3COOH", m = atm$c12$m*2 + atm$h1$m*4 + atm$o16$m*2, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "CH3CN", m = atm$c12$m*2 + atm$h1$m*3 + atm$n14$m, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "CH3OH", m = atm$c12$m + atm$h1$m*4 + atm$o16$m, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "H3PO4", m = atm$h1$m*3 + atm$p31$m + atm$o16$m*4, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "SiO3H2", m = atm$si28$m + atm$h1$m*2 + atm$o16$m*3, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "SiO4H4", m = atm$si28$m + atm$h1$m*4 + atm$o16$m*4, z = 0)),
  list(loss = data.frame(), gain = data.frame(name = "SiC2H6O", m = atm$si28$m + atm$h1$m*6 + atm$o16$m + atm$c12$m*2, z = 0))
)

makeruletable(n1, "n1")



ruleset = rbind(
  makeruletable(isotopes, "iso"),
  makeruletable(z1, "z1"),
  makeruletable(z2, "z2"),
  makeruletable(z3, "z3"),
  makeruletable(n1, "n1")
  )

