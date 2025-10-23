# Contributions
Contributors and maintainers are expected to abide by the [Bioconductor Code of Conduct](https://bioconductor.github.io/bioc_coc_multilingual/en-US.html). 

When creating an issue, use a clear and descriptive title, detail the issue and context in a verbose manner and include code for demonstration. 

We follow a git flow approach for changes. Fork and clone the repository. Create a new branch and develop against current Bioc-devel version. Document documentation changes using devtools::document(). Check and resolve errors and warnings from devtools::check() and BiocCheck::BiocCheck(). Make sure that your branch is synced with the target branch and create a pull request. Detail the changes and context in a verbose manner and tag one of the maintainers for review. Also check the automatic build reports in the pull request.