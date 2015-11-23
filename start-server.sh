#!/bin/bash

R --no-save --vanilla -e "shiny::runApp('.', port=8080, host='0.0.0.0')"
