# FreshmanBioShiny

**Interactive Data Analysis Tool for Understanding Group Behavior Patterns**

**Project Overview:** Developed an interactive R Shiny web application to analyze social behavior patterns from RFID tracking data, enabling users to explore complex network relationships without requiring advanced coding skills. Try the app (https://prashers.shinyapps.io/datascishinyapp/) using the practice datasets found in this GitHub repository. 

**Business Problem Solved:** Traditional data analysis methods for understanding group interactions and space utilization were too complex for end users. This tool democratizes access to sophisticated social network analysis, allowing stakeholders to quickly derive actionable insights from location-tracking data.

**Key Features:**
- Data Integration: Automated merging of RFID detection data with individual identifiers
- Interactive Visualization: Dynamic social network diagrams showing relationship strength and frequency
- Customizable Parameters: Time-window slider allowing users to adjust association criteria based on business needs
- Quantitative Metrics: Automated calculation of network statistics (degree, strength) for each individual
- Export Functionality: CSV download capability for further analysis and reporting

**Technical Implementation:** Built using R Shiny framework with reactive programming for real-time data exploration. Handles temporal datasets and computes network metrics on-demand.

**Impact:** Enables non-technical users to perform sophisticated behavioral analysis, reducing time-to-insight from hours of coding to minutes of interactive exploration. The tool's flexibility makes it applicable to various use cases involving location-based interaction data. It currently **serves 950+ users annually**.

**Skills Demonstrated:** R Programming, Shiny Web Development, Data Visualization, Network Analysis, User Experience Design, Interactive Dashboard Development


**Summary for educators:**
This Shiny app was developed for an introductory biology lab course (BIOL1082L) at the University of Cincinnati to help students analyze and visualize penguin social networks using real RFID detection data. The app supports a newly introduced data science module that familiarizes students with R-based data analysis (app first developed in Fall 2022).

Students participate in interactive lab activities where they simulate penguin movements and generate RFID detection data. They then use the app to explore both their class-generated dataset and real penguin data from the Cincinnati Zoo. The app allows students to:

- Upload RFID detection data and merge it with penguin ID information.
- Adjust time thresholds for defining social associations.
- Generate interactive social network diagrams showing penguin interactions.
- Compute and download social network metrics (degree and strength) for each penguin.

By using this app, students gain hands-on experience with data visualization and network analysis in R without needing extensive coding knowledge.
