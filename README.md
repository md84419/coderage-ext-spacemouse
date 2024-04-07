# CodeRage's coderage.io.spacemouse extension for Omniverse

# Getting Started

## Install Omniverse and some Apps

1. Install *Omniverse Launcher*: [download](https://www.nvidia.com/en-us/omniverse/download)
2. Install and launch one of *Omniverse* apps in the Launcher. For instance: *Code*.

## Option 1: Clone this repo and add this extension to your *Omniverse App*

0. Clone the repo

- `git https://github.com/md84419/coderage-ext-spacemouse.git`

-- or -

- `git clone git@github.com:md84419/coderage-ext-spacemouse.git`

1. In the *Omniverse App* open extension manager: *Window* &rarr; *Extensions*.
2. In the *Extension Manager Window* open a settings page, with a small gear button in the top left bar.
3. In the settings page there is a list of *Extension Search Paths*. Add the `exts` subfolder for this extension there as another search path, e.g.: `C:\projects\coderage-ext-spacemouse\exts`

![Extension Manager Window](/images/add-ext-search-path.png)

4. Now you can find `coderage.io.spacemouse` extension in the top left search bar. Select and enable it.

### A few tips

* Now that `exts` folder was added to the search you can add new extensions to this folder and they will be automatically found by the *App*.
* Look at the *Console* window for warnings and errors. It also has a small button to open the current log file.
* All the same commands work on Linux. Replace `.bat` with `.sh` and `\` with `/`.
* `Extension name` is a folder name in the `exts` folder, in this example: `coderage.io.spacemouse`.
* The most important feature an extension has is a config file: `extension.toml`, take a peek.

## Option 2: Add this extension to your *Omniverse App* without cloning the github repo

Alternatively, a direct link to a git repository can be added to *Omniverse Kit* extension search paths.  Instead of the 'C:\' path above, use this path in the Extension manager (```Extension Manager -> Gear Icon -> Extension Search Path```):

`git://github.com/md84419/coderage-ext-spacemouse.git?branch=main&dir=exts`

Notice that `exts` is the repo subfolder containing the extension(s). More information can be found in ["Git URL as Extension Search Paths"](https://docs.omniverse.nvidia.com/kit/docs/kit-manual/latest/guide/extensions_advanced.html#git-url-paths) section of the [Omniverse developers manual](https://docs.omniverse.nvidia.com/kit/docs/kit-manual/latest/index.html).



## Option 3: Linking with an *Omniverse* app

For a better experience when developing this extension, it is recommended to create a folder link named `app` to the *Omniverse Kit* app installed from *Omniverse Launcher*. A convenience script to use is included.

1. Remove the search path added in the previous section.
2. Open this cloned repo using Visual Studio Code: `code C:\projects\coderage-ext-spacemouse`. It will suggest installing a few extensions to improve the Python experience.

3. In the terminal (CTRL + \`), run:

```bash
> link_app.bat
```

If successful you should see the `app` folder link in the root of this repo.

If multiple Omniverse apps are installed, the script will select the recommended one. Or you can explicitly pass an app:

```bash
> link_app.bat --app create
```

You can also just pass a path to use when creating the link:

```bash
> link_app.bat --path "C:/Users/bob/AppData/Local/ov/pkg/create-2021.3.4"
```

4. Run this app with `exts` folder added as an extensions search path and new extension enabled:

```bash
> app\omni.code.bat --ext-folder exts --enable coderage.io.spacemouse
```

- `--ext-folder [path]` - adds a new folder to the search path
- `--enable [extension]` - enables an extension on startup.

Use `-h` for help:

```bash
> app\omni.code.bat -h
```

5. After the *App* has started you should see:
    * the extension search paths in *Extensions* window as in the previous section.
    * the extension is enabled in the list of extensions.

6. If you look inside `omni.code.bat` or any other *Omniverse App*, they all run *Omniverse Kit* (`kit.exe`). *Omniverse Kit* is the Omniverse Application runtime that powers *Apps* build out of extensions.
Think of it as analagous to `python.exe`. It is a small runtime, that enables all the basics, like settings, Python, logging and searches for extensions. **Everything else is an extension.** You can run only this new extension without running any big *App* like *Code*:


```bash
> app\kit\kit.exe --ext-folder exts --enable coderage.io.spacemouse
```

It starts much faster and will only have extensions enabled that are required for this new extension (look at the `[dependencies]` section of `extension.toml`). You can enable more extensions: try adding `--enable omni.kit.window.extensions` to have the extensions window enabled (yes, extension window is an extension too!):


```bash
> app\kit\kit.exe --ext-folder exts --enable coderage.io.spacemouse --enable omni.kit.window.extensions
```

You should see a menu in the top left. From here you can enable more extensions from the UI.

### A few tips

* In the *Extensions* window, press the *Burger* menu button near the search bar and select *Show Extension Graph*. It will show how the current *App* comes to be: all extensions and dependencies.
* Learn more: [Extensions system documentation](http://omniverse-docs.s3-website-us-east-1.amazonaws.com/kit-sdk/104.0/docs/guide/extensions.html)

# Running Tests

To run tests we run a new process where only the tested extension (and it's dependencies) is enabled. Like in example above + testing system (`omni.kit.test` extension). There are two ways to run extension tests:

1. Run: `app\kit\test_ext.bat coderage.io.spacemouse  --ext-folder exts`

That will run a test process with all tests and then exit. For development mode, pass `--dev`: that will open the test selection window. As everywhere, hot-reload also works in this mode, give it a try by changing some code!

2. Alternatively, in *Extension Manager* (*Window &rarr; Extensions*) find your extension, click on the *TESTS* tab, click *Run Test*

For more information about testing refer to: [testing doc](http://omniverse-docs.s3-website-us-east-1.amazonaws.com/kit-sdk/104.0/docs/guide/ext_testing.html).


# Sharing This Extension to Github

To make the extension available to other users, use [Github Releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository).

1. Make sure the repo has [omniverse-kit-extension](https://github.com/topics/omniverse-kit-extension) topic set for auto discovery.
2. For each new release increment extension version (in `extension.toml`) and update the changelog (in `docs/CHANGELOG.md`). [Semantic versionning](https://semver.org/) must be used to express the severity of API changes.
  - See time index 4:17 in [this video from Mati](https://www.youtube.com/watch?v=lEQ2VmzXMgQ) for more information.

# Known limitations
1. Integration is achieved by simulating (injecting) keyboard and mouse events into the application.
This is a significant limitation as it means that the usual use case of operating the spacemouse with
the non-dominant hand and the regular mouse with the dominant hand, doesn't work - one can only use one
input device at a time.  Hopefully in a future release a better API for injecting Spacemouse commands
will become available.

2. spacenavigator dependency is currently Windows-only. Therefore this extension only works on Windows currently.  There is a fork at https://github.com/JakubAndrysek/PySpaceMouse which we should investigate.
To reduce dependency friction, we may wish to stick with spacenavigator for Windows and only use
pyspacemouse on Linux and Mac OS X

# Contributing
The source code for this repository is provided as-is. We only accept outside contributions from individuals who have signed an Individual Contributor License Agreement.
