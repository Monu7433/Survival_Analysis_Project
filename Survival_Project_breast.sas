/*====================================================
 METABRIC Survival Analysis – SAS
 Author: Monu Prajapati
====================================================*/


/*====================================================
 STEP 1: IMPORT DATA
====================================================*/

options validvarname=v7;

proc import datafile="/home/u64108899/monu/Breast Cancer METABRIC.csv"
    dbms=csv
    out=survival_raw
    replace;

    getnames=yes;
    guessingrows=max;

run;


/*====================================================
 STEP 2: CHECK VARIABLES
====================================================*/

proc contents data=survival_raw;
run;


/*====================================================
 STEP 3: REMOVE MISSING SURVIVAL DATA
====================================================*/

data survival_clean;

set survival_raw;

if missing(Overall_Survival__Months_) then delete;

if missing(Overall_Survival_Status) then delete;

run;


/*====================================================
 STEP 4: CREATE TIME AND STATUS VARIABLES
====================================================*/

data survival_final;

set survival_clean;

time = Overall_Survival__Months_;

if strip(Overall_Survival_Status)="Deceased" then status=1;
else if strip(Overall_Survival_Status)="Living" then status=0;

run;


/*====================================================
 STEP 5: VERIFY DATA
====================================================*/

proc freq data=survival_final;

tables status;

run;

proc means data=survival_final;

var time status;

run;


/*====================================================
 STEP 6: BAR CHART
====================================================*/

proc sgplot data=survival_final;

vbar Overall_Survival_Status;

title "Bar Chart of Survival Status";

run;


/*====================================================
 STEP 7: KAPLAN-MEIER OVERALL
====================================================*/

proc lifetest data=survival_final
plots=(survival(cl atrisk));

time time*status(0);

title "Kaplan-Meier Overall Survival";

run;


/*====================================================
 STEP 8: KM + LOG RANK TEST (ER STATUS)
====================================================*/

proc lifetest data=survival_final
plots=(survival(cl atrisk));

time time*status(0);

strata ER_Status;

title "KM and Log-rank by ER Status";

run;


/*====================================================
 STEP 9: KM + LOG RANK TEST (CHEMOTHERAPY)
====================================================*/

proc lifetest data=survival_final
plots=(survival(cl atrisk));

time time*status(0);

strata Chemotherapy;

title "KM and Log-rank by Chemotherapy";

run;


/*====================================================
 STEP 10: MULTIVARIABLE COX REGRESSION
====================================================*/

proc phreg data=survival_final;

class ER_Status(ref="Negative")
      Chemotherapy(ref="No")
      HER2_Status(ref="Negative");

model time*status(0)=

      Age_at_Diagnosis
      Tumor_Size
      ER_Status
      Chemotherapy
      HER2_Status;

title "Multivariable Cox Regression";

run;


/*====================================================
 STEP 11: PH ASSUMPTION TEST
====================================================*/

proc phreg data=survival_final;

class ER_Status(ref="Negative")
      Chemotherapy(ref="No")
      HER2_Status(ref="Negative");

model time*status(0)=

      Age_at_Diagnosis
      Tumor_Size
      ER_Status
      Chemotherapy
      HER2_Status;

assess ph / resample;

title "Proportional Hazards Assumption Test";

run;


/*====================================================
 STEP 12: STRATIFIED COX MODEL
====================================================*/

proc phreg data=survival_final;

class Chemotherapy(ref="No");

model time*status(0)=

      Age_at_Diagnosis
      Tumor_Size
      ;

strata Chemotherapy;

title "Stratified Cox Model";

run;

proc phreg data=survival_final;

class Chemotherapy(ref="No");

model time*status(0)=
      Age_at_Diagnosis
      Tumor_Size;

strata Chemotherapy;

assess ph / resample;

title "PH Assumption Test - Stratified Cox Model";

run;


/*====================================================
 END
====================================================*/