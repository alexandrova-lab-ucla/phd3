#!/usr/bin/env python3

import logging

from phd3 import scan_coordinates
from phd3 import scan_to_xyz

def scan_2_xyz():
    scan_to_xyz()

def main():
    scan_coordinates()


if __name__ == "__main__":
    main()
