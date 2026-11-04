#!/bin/bash
# Nightly sync of this repo with GitHub (cron, 04:10): commit new/changed files, pull, push, all over
# origin (SSH: git@github.com:sk999cam/ASME_2027.git, no credential prompt). Overleaf pulls from and
# pushes to the same GitHub repo, so the pull also brings in edits made on Overleaf. Each step is logged
# with its result in git_sync.log; the first failure stops the run. Same scheme as 4dns_filmcooled.
REPO="/home/sk999/Documents/ASME_2027"
LOGFILE="$REPO/git_sync.log"
cd "$REPO" || exit 1
export GIT_SSH_COMMAND="ssh -o BatchMode=yes -o ConnectTimeout=30"   # never wait for a prompt

log() { echo "$*" >> "$LOGFILE"; }
step() {   # step <name> <command...>: run, log output and result, stop on failure
    local name=$1; shift
    if timeout 600 "$@" >> "$LOGFILE" 2>&1; then
        log "[ok] $name"
    else
        log "[FAILED] $name (exit $?); sync stopped"
        exit 1
    fi
}

log "=== $(date) ==="
[ "$(git rev-parse --abbrev-ref HEAD)" = main ] || { log "[FAILED] not on main; sync stopped"; exit 1; }

# Order: commit local changes, then pull (merge), then push. Committing first keeps local edits safe:
# git refuses to pull over uncommitted edits to a file that also changed on GitHub (e.g. via Overleaf).

# .gitignore decides what is tracked (template/, LaTeX build files, raw data are excluded)
git add -A >> "$LOGFILE" 2>&1
# GitHub refuses files over 100 MB and Overleaf limits uploads to 50 MB: leave anything over 50 MB unstaged
git diff --cached --name-only --diff-filter=AM -z | while IFS= read -r -d '' f; do
    if [ "$(stat -c %s "$f")" -gt 50000000 ]; then
        git reset -q -- "$f"
        log "[skipped] $f is over 50 MB, not committed"
    fi
done
if git diff --cached --quiet; then
    log "nothing to commit"
else
    step commit git commit -q -m "auto-sync $(date '+%Y-%m-%d %H:%M')"
fi

# a merge, not a rebase: never rewrites commits; a conflict is aborted (local commit kept) and reported
if git ls-remote --exit-code --heads origin main > /dev/null 2>&1; then
    if ! timeout 600 git pull --no-rebase --no-edit origin main >> "$LOGFILE" 2>&1; then
        git merge --abort >> "$LOGFILE" 2>&1
        log "[FAILED] pull (merge conflict with GitHub/Overleaf, merge aborted, local commit kept); sync stopped"
        exit 1
    fi
    log "[ok] pull"
else
    log "origin has no main branch yet: pull skipped"
fi

# also pushes the commits made during the day
step push git push -u origin main
step fetch git fetch origin
log "in sync: $(git rev-parse --short main) = origin/main $(git rev-parse --short origin/main)"
