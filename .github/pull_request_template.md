## Summary

<!-- What this pull request changes and why. -->

## Checklist

- [ ] The base of this pull request is `develop`. Only a promotion from `develop` or an urgent fix (`hotfix/`) targets `main`.
- [ ] The branch was created from the latest base branch and contains its latest commits, so the branch and the base run the same CI.
- [ ] The pull request contains a single, delimited piece of work.
- [ ] The checks pass both on the push of the branch and on the pull request (GitHub Actions and CircleCI), with no failure, error
      or critical skip.
- [ ] Scientific results are unchanged, or every change is described and was approved before being committed.
- [ ] Changes to the method, figures, documentation or interfaces were inspected and approved before being committed.
- [ ] Nothing is published: no tag, release, PyPI, TestPyPI or GitHub Pages step is part of this pull request.

For a promotion from `develop` to `main`, also record the scope, the version, the tests, the regressions and the scientific impact.
A tag is created only after the version and the tag are explicitly authorized.
