# SESSION 1: HPCC INITIATIVE


---

Welcome to the Texas Tech University 2021 Bioinformatics Workshop, held by the TTU Chapter Botanical Society Association. The first session of the series will guide you through the very first steps of using Texas Tech's High Performance Computer Cluster and navigating yourn personal space within the cluster. We will finish off by downloading sequence data from the KEW database.



---

[TOC]

## 1.1: Opening Your Terminal

All computers allow you access to a terminal, which is a window that allows you to communicate with your computer, and command it to perform functions. All computers have a built in terminal. 

In Windows, you can click the Start Icon in the bottom left on your screen. In the search field, you can type 'terminal', 'cmd', or 'Command Prompt' and the application Command Prompt will appear. Alternatively, you can access Windows Powershell by right clicking on the Start Icon, and finding the application on the menu. Powershell is Windows' "modern command shell that includes the best features of other popular shells".

On a Mac, there are two methods to opening your terminal. First, you can use the Finder by clicking on the Finder logo, clicking Applications, opening the Utilities folder, and double clicking "terminal". You could also find the terminal by pressing the Command buttom and the space bar simultaneously. Type 'terminal' into the search bar, and then double click the terminal application to launch it.

Many scientists who frequently use command line to process data for their research utilize "Gitbash", which can be downloaded [here](https://git-scm.com/downloads). Gitbash is available for macOS, Windows, and Linux/Unix system. The advantages that gitbash offers over your personal computer's terminal is that it understands a wider range of commands, making certain projects easier. 


## 1.2 Logging into the Red Raider Cluster

You will be utilizing the new ['RedRaider' cluster](https://www.depts.ttu.edu/hpcc/about/Introduction-to-RR-January-2021.pdf), which is a new, high powered cluster available to students at TTU. To log onto the cluster, you must be considered 'on-campus'. This can be achieved two ways;

1. Use the wired TTU network or the TTUnet wireless network. Connections from TTUHSC or the TTUguest network are not considered 'on-campus'.
2. If you are physically off-campus, use the [TTU Global Protect VPN client](https://www.askit.ttu.edu/portal/app/portlets/results/viewsolution.jsp?&solutionid=140702103827226).

*Off-campus connection is possible, but less convenient. For a walkthrough from accessing the cluster while 'off-campus', follow [this guide](https://www.depts.ttu.edu/hpcc/userguides/generalguides/logingeneral.php#linux-offcampus) which can be found on the TTU HPCC website. If you are having trouble with 'on-campus', [this page](https://www.depts.ttu.edu/hpcc/userguides/general_guides/login_general.php#linux-oncampus) also explains logging on in greater detail*

Once you have satisfied 'on-campus' requirements, you will type this command into the command line:

```
ssh <eraider>@login.hpcc.ttu.edu
```
Hit the enter key to run the line. You will then be prompted to input your eraider password. Once initialized, you will be connected to the cluster and able to run commands. 
## 1.3 Working in an Interactive System

The 'interactive' command allows you to directly access a CPU or GPU. This allows the user to allocate resources. Nodes can also be directly interacted with. Additionally, by using an interactive session, results from our commands are displayed on the screen.


HPCC has outlined the command and arguments below:
> Command: 
> interactive [-A] [-c] [-p] [-J] [-w] [-g] [-h]
> 
> Optional arguments:
> 
> -A: the account name
> 
> -c: number of CPU cores to request (default: 1)
> 
> -p: partition to run job in (default: nocona)
> 
> -J: job name (default: INTERACTIVE)
> 
> -w: node name
> 
> -g: number of GPU to request
> 
> -h: show this usage info

For, now, lets use the default command for this portion of the workshop.

```
interactive
```

## 1.4 Making a Hierarchy File System
Lets try a few commands that allow you to visualize and organize your personal directory in command line.

> ls   :     lists all directories in the directory you are currently in. Use 'ls' to identify folders in your current directory

> mkdir     : Makes a new directory in the current directory you are in. 
> 

Go ahead and create a new directory to keep your data in for this tutorial. You can call it whatever you'd like, for tutorial purposes, we will call our's TC_Workshop. After youve created it, use the ls command to see it in your files.

`mkdir [newname]`

> cd  : Change Directory. You can use this command to move between directories in your terminal.

Navigate to your new directory by using the cp command.
```
cd [newname]
```
You can return to your previous directory by using the entire path to the directory, of by using '..' which denotes the directory that yout current directory is in. Try moving around your terminal using these commands.

```
cd ../
OR
cd /scratch/[username]
```
```
cd /TC_Workshop
```
If you get lost, try using the list command (ls) to get a hint about where you can move.
## 1.5 Copy Data Across Folders

To demonstrate how to copy data between two folders, we will create a text file and copy it into an adjacent directory. Return to your directory that you just created (TC_Workshop). For the purposes of this tutorial, we will call our directory TC_Workshop. Make two new directories within the TCWorkshop directories, called Test1 and Test2.

```
mkdir Test1
mkdir Test2
```
Navigate into the Test1 directory.
```
cd /Test1
```
To create a new text file, you can use the nano command.
```
nano
```
This command opens the default nano screen. Type any sort of message into the screen. To exit nano, use Control X. This then prompts you to name the file. If you use the 'ls' command, you can now see a text file in your Test1 directory. 

Navigate to your Test2 directory using cd commands. Once you are in Test2, you can utilize the cp command to copy your textfile into this folder.

```
cp /scratch/username/TC_Workshop/Test1/filename.txt /scratch/username/TC_Workshop/Test2
```
to simplify this command further, we can use '.'. We used the double period ".." to define the parent directory earlier. Try using "." instead of the filepath to the Test2 directory.

If you use the ls command in the Test2 directory, you should see a copy of the text file. You can open it with nano to see a copy of the same message you wrote previously. The cp command can be used across directories in your space, and other user's space on HPCC, as long as you have the pathname to retrieve the file.
## 1.6 Opening and Modifying Your /bashrc File

## 1.7 Running Singularity Containers

## 1.8 Downloading data from KEW database
