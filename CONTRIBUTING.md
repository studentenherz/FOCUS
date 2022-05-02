# Contributing to FOCUS

First of all, thanks for taking the time to contribute! 🎉

Here are some guidelines for contributing to FOCUS, they are more suggestions than rules, so use your best judgement and feel free to propose changes to this document.

### Keeping the code portable

The following guidelines aim at keeping the codebase clean and easy to understand, free of any snippet of code or files that are specific to the local environment in which the code is being used and modified, or the case being studied.

- **No compiled objects** should be in the repository, only code.
- Input data should be detached from code. **Input data should only exist** in the repo **for testing purposes** inside the `test` directory.
- **No local configuration files** (i.e.: IDE generated config, build folders, etc.) should be in the repo, make sure to populate your `.gitignore` in order to exclude them.
- There should not be any absolute path inside the code.
- Required system properties should be obtained at runtime or during compilation, never hard-coded.
- Any data required to _run an experiment_ should be obtained at runtime through configuration files, input files, etc., never hard-coded into the source.
- The **user shouldn't need to change the code** except to make improvements or fix bugs.

### Git Commit Messages

When committing code to the project provide a good message explaining the change. In order to create good commit messages follow **the seven rules** (deeper explanation [here](https://cbea.ms/git-commit/)):

1. Separate subject from body with a blank line.
2. Limit subject line to 50 characters.
3. Capitalize the subject line.
4. Do not end the subject line with a period.
5. Use the imperative mood in the subject line.
6. Wrap the body at 72 characters.
7. Use the body to explain _what_ and _why_ instead of _how_.

### Coding Style

- Use descriptive names for variables and functions. Some seconds thinking of a good name could save minutes of code reading later.
- **Use English**, both for the comments and the code.
- Prefer tabs over spaces for indentation purposes. It saves disk space.
- In C/C++, provide closing braces with a comment when the block is large, for example:
  ```C++
  namespace focus {
  	// many many lines of code
  } // namespace focus
  ```
- Document the functions with Doxygen compatible comments in the format of the following example:

  ```C++
  /* Adds two doubles and returns the result.
   * @param x first number.
   * @param y second number.
   * @return the sum of the two numbers.
   */
  double sum(double x, double y){
  	return x + y;
  }
  ```

  See more [here](https://www.doxygen.nl/manual/docblocks.html). Some editors (like vscode, for instance) can show this documentation as tooltips while writing code.

Can't think of anything else right now, just keep your code clean and make useful comments. Maybe refer to the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
