#!/bin/bash
for N in 252 250 253 254; do echo "=== $N dense ==="; ssh $N 'tail -2 ~/dense_run/dense.log'; done
for N in 255 251; do echo "=== $N magnetic ==="; ssh $N 'tail -2 ~/mag_run/mag.log'; done
