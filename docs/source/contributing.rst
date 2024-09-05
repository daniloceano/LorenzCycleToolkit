Contributing
============

We welcome contributions to improve LorenzCycleToolkit! Please follow these steps to contribute:

1. **Fork the repository** and clone your fork.
2. **Create a new branch**: 
   - Use one of the following naming conventions:
     - ``feature/your-feature-name`` for new features
     - ``bugfix/your-bug-name`` for bug fixes
     - ``hotfix/your-hotfix-name`` for urgent fixes
   - Example: ``git checkout -b feature/new-visualization``
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
      
8. **Create a pull request**: Submit your PR through GitHub, following our PR template and referencing any related issues.

Before submitting your pull request, please ensure:

- Your code passes all tests and follows the coding standards.
- You have added or updated documentation if necessary.

Thank you for contributing!
