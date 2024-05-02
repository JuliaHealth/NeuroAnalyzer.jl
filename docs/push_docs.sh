#!/bin/bash

GIT_STATUS=`git status | awk 'FNR==2'`
if [ "$GIT_STATUS" != "Your branch is up to date with 'origin/main'." ]; then
  git commit -m "Woodpecker CI Documentation Build"
  git push -f origin main
fi