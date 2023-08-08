import os
import pandas
import sys
import yaml
import urllib.request

shell.executable("bash")

# default configuration file.
configfile:
    srcdir("config/config.default.yaml")

# some parameters.
