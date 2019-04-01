#!/bin/bash
echo $HOST > host.txt

vasp-ver-bsub-native bfscan -n 16 -o qlog.txt -q suncat3 -W 20:00
