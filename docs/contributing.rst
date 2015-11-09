============
Contributing
============

Contributing knowledge
----------------------

If you've found a bug, or have a suggestion for improving the toolkit,
please let us know by creating a GitHub issue at:

https://github.com/ParkerLab/atactk/issues


Contributing code
-----------------

We welcome contributions of code, too. Here's how to get started.

#. Fork the repo on GitHub: https://github.com/ParkerLab/atactk/

#. Set up a Python virtualenv. Assuming you have virtualenvwrapper
   installed, this is how you set up your fork for local development::

     $ mkvirtualenv atactk
     $ git clone git@github.com:your_name_here/atactk.git
     $ cd atactk/
     $ python setup.py develop

#. Create a branch for local development::

     $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

#. `(Optional, but much appreciated.)` When you're done making changes,
   check that your changes pass flake8 and the tests, including
   testing other Python versions with tox::

     $ flake8 atactk tests
     $ python setup.py test
     $ tox

   To get flake8 and tox, just `pip install` them into your virtualenv.

#. Commit your changes and push your branch to GitHub::

     $ git add .
     $ git commit -m "Your detailed description of your changes."
     $ git push origin name-of-your-bugfix-or-feature

#. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

#. If the pull request adds functionality, please update the
   documentation, especially docstrings and script help text.

#. If you have time to write tests too, great, but we understand
   you're volunteering your time to help our project, and we will
   take care of making sure changes are tested.
