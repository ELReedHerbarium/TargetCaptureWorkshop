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


## 1.4 Making a Hierarchy File System

## 1.5 Copy Data Across Folders

## 1.6 Opening and Modifying Your /bashrc File

## 1.7 Running Singularity Containers

## 1.8 Downloading data from KEW database
