# Installing ProPrep

This guide assumes nothing: no conda, no Python, and no experience with the
command line. Follow the section for your computer. Each one ends with the same
test, so you know it worked.

**What you will download:** one installer file of about 1 to 1.5 GB. It
contains everything ProPrep needs (ProPrep itself, AmberTools, MODELLER and
their dependencies). **What you will need:** about 8 GB of free disk space and
10 to 15 minutes. The install itself needs no internet connection once the
file is downloaded.

**What you will get:** a folder called `ProPrep` in your home directory, and
inside it a program `proprep-web` that opens ProPrep in your web browser.

---

## Step 1: find out which installer you need

| Your computer | Installer file | Section |
|---|---|---|
| Mac with an Apple chip (M1, M2, M3, M4, ...) | `ProPrep-<version>-MacOSX-arm64.sh` | A |
| Mac with an Intel chip | `ProPrep-<version>-MacOSX-x86_64.sh` | A |
| Windows 10 or 11 | `ProPrep-<version>-Linux-x86_64.sh` (runs inside WSL, see below) | B |
| Linux | `ProPrep-<version>-Linux-x86_64.sh` | C |

**Not sure which Mac you have?** Click the Apple menu (top-left corner) and
choose **About This Mac**. If the window says *Chip: Apple M...* you have an
Apple chip; if it says *Processor: Intel...* you have an Intel chip.

Download the file from the ProPrep releases page:

> https://github.com/mjgplab/proprep/releases/latest

Scroll down to **Assets** and click the file name. Your browser will save it
to your **Downloads** folder. The download is large; give it a few minutes.

---

## Section A: macOS (Apple chip or Intel)

### A1. Open the Terminal

The Terminal is the app where you type commands. Press **Command + Space**,
type `Terminal`, and press **Return**. A window with a text prompt appears.
Everything below is typed into this window, one line at a time, pressing
**Return** after each line.

### A2. Run the installer

Type the line that matches your Mac. (Replace `<version>` with the version
number in the file you downloaded, for example `1.17.0`. Tip: type `bash
~/Downloads/ProPrep-` and then press the **Tab** key; the Terminal completes
the file name for you.)

Apple chip:

```
bash ~/Downloads/ProPrep-<version>-MacOSX-arm64.sh -b -p ~/ProPrep
```

Intel chip:

```
bash ~/Downloads/ProPrep-<version>-MacOSX-x86_64.sh -b -p ~/ProPrep
```

The installer prints `Unpacking payload...` and then a list of packages. It
takes 2 to 5 minutes. When it finishes it prints `ProPrep is installed.` and
returns you to the prompt. `-b` means "accept the license and do not ask questions"; `-p ~/ProPrep`
is where it installs: a folder named `ProPrep` in your home folder.

If the Terminal says `Permission denied`, you typed the file name without
`bash` in front. Retype the whole line starting with `bash`.

### A3. Test it

```
~/ProPrep/bin/proprep --version
```

You should see `ProPrep v<version> - Proper Protein Preparation workflow manager`.

### A4. Start ProPrep

```
~/ProPrep/bin/proprep-web
```

Your web browser opens with ProPrep running inside it. Leave the Terminal
window open while you work; closing it stops ProPrep. To stop ProPrep, close
the browser tab (the program shuts down by itself a moment later) or press
**Control + C** in the Terminal.

Now go to **Step 3: MODELLER license key** below.

---

## Section B: Windows 10 or 11

ProPrep runs on Linux and macOS. On Windows it runs inside the **Windows
Subsystem for Linux (WSL)**, a feature built into Windows that runs a small
Linux system alongside Windows. You install WSL once (Section B1), and after
that ProPrep behaves like any other program: it opens in your normal Windows
web browser.

### B1. Install WSL (once)

1. Click the **Start** button, type `PowerShell`, right-click **Windows
   PowerShell** and choose **Run as administrator**. Click **Yes** if Windows
   asks for permission.
2. In the blue window, type the following and press **Enter**:

   ```
   wsl --install
   ```

   Windows downloads and installs WSL and a Linux distribution called Ubuntu.
   This takes several minutes.
3. **Restart your computer** when it asks you to.
4. After the restart, an **Ubuntu** window opens on its own and says
   *Installing, this may take a few minutes...*. It then asks you to choose a
   **username** (lowercase, no spaces) and a **password**. The password does
   not show on screen while you type; that is normal. Remember it.

   If no Ubuntu window opens after the restart, click **Start**, type
   `Ubuntu`, and open it.

You now have an Ubuntu window with a text prompt ending in `$`. This is the
**Ubuntu terminal**; every command below is typed here, one line at a time,
pressing **Enter** after each.

Requirements: Windows 10 version 2004 or later, or Windows 11; about 10 GB of
free disk space in total. If `wsl --install` reports an error, the most
common causes are an older Windows 10 (run Windows Update) or virtualization
disabled in the computer's firmware (BIOS/UEFI); your IT support can enable it.

### B2. Run the installer

Your Windows **Downloads** folder is visible from Ubuntu under
`/mnt/c/Users/<YourWindowsName>/Downloads`. Type this line, replacing
`<YourWindowsName>` with your Windows user name and `<version>` with the
version in the file name (for example `1.17.0`):

```
bash /mnt/c/Users/<YourWindowsName>/Downloads/ProPrep-<version>-Linux-x86_64.sh -b -p ~/ProPrep
```

Tip: type the line up to `Downloads/ProPrep-` and press the **Tab** key; the
terminal completes the file name. If you do not know your Windows user name,
type `ls /mnt/c/Users` and press Enter to see the list.

The installer prints `Unpacking payload...` and then a list of packages. It
takes 2 to 5 minutes. When it finishes it prints `ProPrep is installed.`.
ProPrep is installed inside Ubuntu, in a folder named `ProPrep` in your Ubuntu
home folder (that is what `-p ~/ProPrep` asked for).

### B3. Test it

```
~/ProPrep/bin/proprep --version
```

You should see `ProPrep v<version> - Proper Protein Preparation workflow manager`.

### B4. Start ProPrep

```
~/ProPrep/bin/proprep-web
```

The Ubuntu window prints a line like `ProPrep Web Shell listening on
http://127.0.0.1:8000/`. Open your normal Windows web browser (Edge, Chrome,
Firefox) and go to:

> http://127.0.0.1:8000/

ProPrep appears in the browser. Leave the Ubuntu window open while you work.
To stop ProPrep, close the browser tab (the program shuts down by itself a
moment later) or press **Control + C** in the Ubuntu window.

**Where are my files?** ProPrep saves your projects in the Ubuntu home folder.
To see them in Windows Explorer, type this in the Explorer address bar:

> `\\wsl$\Ubuntu\home\<your-ubuntu-username>`

Now go to **Step 3: MODELLER license key** below.

---

## Section C: Linux

Open a terminal. Then:

```
bash ~/Downloads/ProPrep-<version>-Linux-x86_64.sh -b -p ~/ProPrep
~/ProPrep/bin/proprep --version
~/ProPrep/bin/proprep-web
```

The installer takes 2 to 5 minutes and puts ProPrep in `~/ProPrep`
(`-b` skips the questions, `-p` sets the location).
`proprep-web` opens your default browser; on a machine without a graphical
desktop it prints the address to open instead. The installer is built for
x86-64 Linux and has been tested on Rocky Linux.

Now go to **Step 3: MODELLER license key** below.

---

## Step 3: MODELLER license key (needed for structure repair)

ProPrep uses MODELLER to fill in missing parts of protein structures. MODELLER
is free for academic use but requires a personal license key.

1. Register at https://salilab.org/modeller/registration.html (name, email,
   institution; choose the academic license). The key arrives by email
   within a few minutes. It looks like a short string of letters and digits.
2. In your terminal (Terminal on a Mac, the Ubuntu window on Windows), type
   the following, replacing `YOUR_KEY` with the key from the email, keeping
   the quotation marks:

   ```
   mkdir -p ~/.proprep && echo 'YOUR_KEY' > ~/.proprep/modeller_key
   ```

That is all. ProPrep reads the key every time it starts. You never need to
enter it again, and it is never sent anywhere. If you mistype the key, run the
same line again with the correct one.

Without a key ProPrep still runs; only the structure-repair step is disabled,
and ProPrep says so when it starts.

---

## Checking that everything works

Run these two lines (Terminal on a Mac, Ubuntu window on Windows):

```
~/ProPrep/bin/proprep --version
~/ProPrep/bin/proprep-web
```

The first prints the version. The second opens ProPrep in your browser (on
Windows, open http://127.0.0.1:8000/ yourself). If ProPrep's start-up text
does not mention MODELLER being unavailable, your license key was found.
Press **Control + C** in the terminal to stop.

---

## Problems

**`No such file or directory`** when running the installer: the file name or
folder is wrong. Use the Tab key to complete the name rather than typing it
out, and check that the download actually finished (the file should be about
1 GB or more).

**`Permission denied`**: the line must begin with `bash`.

**`ERROR: File or directory already exists: ~/ProPrep`**: ProPrep is already
installed. To reinstall or update, add `-u` to the command (see "Updating"
below).

**The installer stops with `No space left on device`**: free up at least 8 GB
and run it again.

**On a Mac, macOS asks whether to allow Terminal to access your Downloads
folder**: click **OK**; the installer only reads the file you downloaded.

**On Windows, `wsl --install` says it needs a restart, then nothing happens**:
open **Start**, type `Ubuntu`, and run it; the first-time setup continues from
there.

**The browser shows "This site can't be reached" at 127.0.0.1:8000**: ProPrep
is not running. Start `~/ProPrep/bin/proprep-web` in the terminal and try
again. If the terminal says port 8000 is busy, it will pick 8001; use the
address it prints.

**Antivirus or corporate software blocks the installer**: the `.sh` file is a
plain script that unpacks files into your home folder; ask your IT support to
allow it, or install on a personal computer.

## Removing ProPrep

Delete the `ProPrep` folder in your home directory (`rm -rf ~/ProPrep` in
the terminal). Your projects and your MODELLER key are stored elsewhere and
are not removed.

## Updating to a newer release

Download the newer installer from the releases page and run it as above with
one extra option, `-u`, which tells it to update the existing folder:

```
bash ~/Downloads/ProPrep-<version>-<platform>.sh -b -u -p ~/ProPrep
```

Your projects and your MODELLER key are stored elsewhere and are untouched.
(Without `-u` the installer stops with `ERROR: File or directory already
exists`.)
