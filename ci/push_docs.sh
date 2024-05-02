#!/bin/bash

GIT_STATUS=`git status | awk 'FNR==2'`
[ "$GIT_STATUS" != "Your branch is up to date with 'origin/main'." ] && git push -f origin main
