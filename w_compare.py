#!/usr/bin/env python

import numpy as np
import cPickle as pickle

#w0 = pickle.load(open('VAup23-26_lin.pk', 'rb'))*3600.0
#w1 = pickle.load(open('VAup24-27_lin.pk', 'rb'))*3600.0
w0 = pickle.load(open('tendxall23-26_lin.pk', 'rb'))
w1 = pickle.load(open('tendxall24-27_lin.pk', 'rb'))

print w0[1,:,110,:].max(), w1[0,:,110,:].max()
