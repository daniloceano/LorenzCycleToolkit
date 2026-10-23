.. EDITORIAL PROPOSAL (2026-10-23), awaiting the project lead's review: steps 2 and 8, the last item of the list before
   submitting and the section "Branches and releases" follow the branch policy in docs/development_workflow.md.

Contributing
============

We welcome contributions to improve LorenzCycleToolkit! Please follow these steps to contribute:

1. **Fork the repository** and clone your fork.

2. **Create a new branch** from the latest ``develop`` branch of this repository (a ``hotfix/`` branch starts from ``main``):

   - Use one of the following naming conventions:
   
     - ``feature/your-feature-name`` for new features
     - ``bugfix/your-bug-name`` for bug fixes
     - ``hotfix/your-hotfix-name`` for urgent fixes
   
   - Example, with the remote ``upstream`` pointing to ``daniloceano/LorenzCycleToolkit``:

     .. code-block:: bash

        git fetch upstream
        git checkout -b feature/new-visualization upstream/develop

3. **Make your changes** following the project's coding standards.

4. **Test your changes**: Ensure all tests pass by running:

   .. code-block:: bash

      pytest

5. **Format your code**: Run the following commands to ensure your code adheres to the project standards:

   .. code-block:: bash

      autopep8 --in-place --recursive .
      isort .
      flake8

6. **Commit your changes**: 
   - Use descriptive commit messages, e.g., ``git commit -m 'Add energy cycle visualization feature'``
7. **Push your changes** to your fork:
   .. code-block:: bash

      git push origin feature/new-visualization
      
8. **Create a pull request** against the ``develop`` branch (a ``hotfix/`` pull request targets ``main``): Submit your PR through GitHub, following our PR template and referencing any related issues.

Before submitting your pull request, please ensure:

- Your code passes all tests and follows the coding standards.
- You have added or updated documentation if necessary.
- Your branch contains the latest commits of its base branch.

Branches and releases
---------------------

Pull requests are merged into ``develop``. The maintainers promote ``develop`` to ``main`` by a separate pull request when a
set of changes is ready, and an urgent fix merged into ``main`` is then merged into ``develop``. Tags are created on ``main``,
and publishing to PyPI, TestPyPI or GitHub Pages is a separate, manual step. The full policy is in
`docs/development_workflow.md <https://github.com/daniloceano/LorenzCycleToolkit/blob/develop/docs/development_workflow.md>`_.

Thank you for contributing!
