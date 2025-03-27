# CG_A2_202112354

This repository contains a solution with three projects. Follow the instructions below to download the solution, set up the startup project, and run all three projects.

## Prerequisites

Before you begin, ensure that you have the following software installed on your system:

- [Git](https://git-scm.com/) (for cloning the repository)
- [Visual Studio](https://visualstudio.microsoft.com/) (for building and running the projects, 2022 is required)
- [OpenGL](https://www.opengl.org/) (graphics rendering library)
- [GLEW](http://glew.sourceforge.net/) (OpenGL extension management library)
- [GLFW](https://www.glfw.org/) (window creation and input management library)
- [GLM](https://github.com/g-truc/glm) (OpenGL mathematics library)

This project requires the same additional libraries as the provided sample code for the assignment.


## Cloning the Repository

To get started, clone the repository to your local machine using Git Bash.

1. Open Git Bash on your computer.
2. Navigate to the directory where you want to clone the repository.
3. Run the following command to clone the repository:
    git clone https://github.com/asbbicsar/CG_A2_202112354.git
4. Once the repository is cloned, navigate into the project folder:
    cd CG_A2_202112354


## Opening the Solution in Visual Studio

1. Open the solution file(CG_A2_202112354.sln) in Visual Studio.
   You can do this by double-clicking the .sln file, or by opening Visual Studio and selecting "Open a Project or Solution" and navigating to the solution file.


## Change the Startup Project
To run all the projects, you need to set one of them as the startup project.

1. Right-click on the project you want to set as the startup project in Solution Explorer.
2. Choose Set as StartUp Project.
    If you want to run multiple projects simultaneously, you can configure multiple startup projects:
    - Right-click on the solution in Solution Explorer and select Properties.
    - In the Common Properties section, select Startup Project.
    - Choose Multiple Startup Projects and set the action for the projects you want to run to Start.
    - This way, all selected projects will start.


## Running the Projects
After setting the startup project, you can run the projects using the following method:

1. Press Ctrl+F5 key or click 'Start Without Debugging' in Visual Studio to run the startup project.
2. If you have configured multiple startup projects, all selected projects will run simultaneously.

---
This repository is developed with assistance from ChatGPT.
