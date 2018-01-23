# Incubating-a-DREAM
This is a paper-project repository exploring the optimal use of lapse rate physics in earth surface models supported by the University of Washington eScience Institute Winter Incubator project (2018) and the National Science Foundation PREEVENTS project led by the UW Civil & Environmental Engineering Department.

Visit
1. [Incubating-a-DREAM wiki](https://github.com/ChristinaB/Incubating-a-DREAM/wiki) to see the paper in progress. 
2. Prezi Presentation for [Visuals and Figures](https://prezi.com/view/X9KGi7p9zhEchXHwugcH/).
3. HydroShare Resources
  * [DHSVM Model Instance for Sauk Watershed](https://www.hydroshare.org/resource/abe2fd1e3fc74e889b78c2701872bd58/)
  * [DHSVM Forcing Files for Sauk Watershed](https://www.hydroshare.org/resource/9db8e06bf0024dc9a8aaabbcecc37416/)
  * [Observed Streamflow Download for Sauk River](https://www.hydroshare.org/resource/1ddfc8fd736f4947a04b8e0620991854/)
  * [Jupyter Notebooks for DHSVM + Landlab + DREAM Incubator Workflows](https://www.hydroshare.org/resource/7c3416535ab24d4f93b0b94741bb9572/)

# Vision 
Sparcity of temporal data in high elevations is a hurdle for every water resources researcher investigating hydrologic processes in watersheds with greater than 1000 m rise in elevation.  Temperature is interpolated to high elevations from low elevation measurements as a lapse rate (C/km), are correlated to precipitation (rain or snow) and short and long wave radiation.  Across the continental United States, an annual average temperature lapse rate of 6.5 C/km is commonly used; in the Pacific Northwest, where clear-sky assumptions do not hold, an annual average temperature lapse rate of 4.8 C/km has been reported (Minder et al, 2012).  Although annual average constant lapse rates are representative for the long term average hydrologic response, the sensitivity and importance of sub-daily, sub-watershed distributed mountain microclimatology and atmospheric mechanisms that control temperature lapse rates for capturing atmospheric river, rain on snow, and snowmelt driven flood events is not well understood. Interpolated and gridded shorter time scale lapse rates (e. g., daily) based on low elevation weather observations (COOP stations) are generated at continental scales with no local control on lapse rate assumptions. We need to understand the tradeoff between capturing long term watershed behavior (annual or monthly average) and short term events dependent on high elevation atmospheric interactions (with limited observations).

To address this problem we plan to develop a process to update gridded hydrometerology climate forcings at a watershed scale using local datasets, interpolating across elevation ranges based on theoretical and empirical relationships (e.g. Unsworth and Monteith, 1975). We will use data collected beginning with 2015 by the researchers from the UW Watershed Dynamics Research Group and the Nooksack Indian Tribe. The temperature observation field campaign is across a North Fork Nooksack transect of Mt. Baker (600-1800 elevation range). These results show that 3hr time series of lapse rates vary +/- 5 C/km from an 4.8 C/km annual average based on precipitation (or clouds indicated by relative humidity. Our proposed methods are to use Dakota software (Adams et al, 2015), and a multi-objective optimization algorithm (DREAM; Vrugt, 2016) to calibrate a Landlab (Hobley et al, 2017) Hydrometeorology component that uses low elevation observations (following Livneh et al, 2015) and physics-based lapse rates to match high elevation atmospheric model (WRF; Henn et al., 2017; Currier, 2016; Currier et al., 2017) monthly averages, from a cloud computing environment (e.g. HydroShare).  The feasibility of our proposed process is supported by the fact that our team and projects (PREEVENTS, Landlab, HydroShare) are published authors in each field of expertise; during this incubator we will combine existing utilities to address this problem. 

# Objectives
In our improved process, any researcher with Python coding skills will be able to efficiently generate updating hydrometeorology forcings using local datasets, locally validated assumptions, anywhere in the world, for any time period of interest. We propose a 4-step research plan:

Aim 1 will be to conduct the data wrangling for observed and modeled hydrometeorology time series. Long term climate stations and WRF gridded hydrometeology (Salathe et al., 2014) will be used to assess the lapse rates between the low elevation observations and the long term physics based monthly averages at high elevations.

Aim 2 will be to develop a Landlab component that can flexibly reproduce a grid of interpolated hydrometeorology variables. This will include using empirical and theoretical relationships to derive high elevation correlations between precipitation, temperature and radiation using observations at lower elevations.

Aim 3 will use DREAM multiobjective optimization in Dakota to calibrate empirical parameters.  This will generate a gridded time series of hydromet variables that match WRF long term averages, updated to near real time, and with correlations that can be used to generate long time series (1915-2099) based on corrections to existing gridded datasets (Livneh et al., 2013;  MACA: Abatzoglou and Brown, 2012). 

Aim 4 will be to execute the process in a cloud environment.  The Consortium of Universities Allied for Hydrologic Science (CUAHSI) has a docker image of the software environment needed for us to use the proposed tools. We propose to work with UW Cloud Computing to launch the process from HydroShare, using a dockerized image on a commercial cloud platform to compare the costs and benefits between using CyberGIS (currently used), Azure and Amazon Web Services.

This project establishes a process to overcome current limitations in continental and global scale hydrometeorology datasets to provide greater control over the use of local datasets and microclimatology to improve interpolation of hydrometeorology across data sparse mountain watersheds. The open source Landlab component and use of multi-objective optimization will provide a valuable resource for the wider scientific community struggling to select which climate forcing dataset is best for representing local watershed scale processes - from high mountain Asia, to South America cordillera, and to the North Cascades in the Pacific Northwest.	

# Success Criteria

Reduce Q residuals using DREAM
Test lapse rate sensitivity and optimization

Completed experimental design for a target audience and focused study

# Deliverable Timeline

Week 1-2 Refine Narrative and Milestone Diagram

Week 2-4 Build Landlab Component, Run DREAM for existing DHSVM in Skagit watershed

Week 5-6 Refine code and strategize around obstacles

Week 7-8  DREAM with more experiments

Week 9-10 Draft Paper, Prepare Code and Data for Publication
