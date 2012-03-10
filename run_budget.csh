#!/bin/csh

foreach s (10 16 22 28 34 40 46 52 58 64)
	@ e = $s + 5
	echo "./mom_budget.py $s $e"
	./mom_budget.py $s $e
end
