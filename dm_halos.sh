#!/bin/bash

python3 undimensionalize_changeunits.py
python3 g_swave.py
python3 g_pwave.py
python3 g_dwave.py
python3 g_som.py
python3 plotting.py
