#Data cleaning

General and simplified flow of analysis
1. Purpose. The purpose must be well-defined to let the analysts know what information should they extract from the data. The purpose will affect the type of model to be used and therefore the data preparation process. For example, if prediction accuracy is the ultimate concern of the stakeholders instead of trying to understand the underlying relationships, complex models, such as boosted regression tree, random forest, and SVM can be adopted. These models may give high prediction accuracy, but difficult to understand the underlying relationships between variables. 
2. Data preparation. 
3. Modelling. 
4. Reporting. 

Organising ICD 9 and 10 code with stringr

ICD is the abbreviation of International code of diseases. It is commonly used to document diseases with a set of numbers (ICD 9) and an alphabet with numbers (ICD 10). It is used in calculating healthcare expenditure, demand, and co-morbidity scoring. 

The purpose of this article is to group both ICD 9 and 10 codes into different groups with the help of stringr. The groups of cardiovascular diseases are described in Lagerweij et al., 2017. 

ICD for identifying CVD subgroups in FRS 
Coronary Heart disease -> IHD ( I20-I25 )
Stroke -> Cerebrovascular disease (I60-I69)
Congestive heart failure -> Heart failure (I50)
Intermitten Claudiation -> Atherosclerotic Arteries disease (I70)

Code of concern (): 
MI: I21, I22 ; 
OCHD: I20, I23, I24, I25; 
Cardiac arrest: I46, R96; 
Hemorrhagic stroke: I60, I61, I62; 
Ischemic stroke: I63, I65; 
Other Stroke: I64, I66; 
Other CVD: G45, I67, I69, I70-I74, I50; 
To add: TIA (G45), R96

pacman::p_load(dplyr, lubridate, tidyr, stringr,Hmisc) #load all package required

hos_r <- spss.get("../Data/hos_CVD.sav") 

Hospitalisation data in spss file format is loaded into R. 

hos_raw <- hos_r %>% 
  mutate(date_h = ymd(INDATUMA)) %>% 
  subset(select=c(LopNr, date_h, HDIA3)) %>%
  arrange(LopNr) %>%
  distinct(LopNr, .keep_all = TRUE) #remove duplicated entries
  
Few things has happened with the above code. 
1. A new dataset is created from the imported dataset. 
2. The date is tidied with lubridate and put into the format year-month-day. 
3. There were many different variables in the dataset, but only 3 are important to the study at this moment. So only 3 columns are kept. 
4. LopNr is an identity code. The dataset is arranged based on the rank of LopNr. It is for the easy of later comparison, since individual could have entered into the hospital more than once. 
5. Since we only need the first every hospitalisation, only the first is kept with the distinct function. The .keep_all is to keep all other columns after removing the duplicates. 

hos_icd <- hos_raw %>% 
  mutate(icd = str_trim(as.character(droplevels.factor(HDIA3)), "right")) %>% 
  subset(select = c(LopNr, date_h, icd)) %>%
  mutate(group_h = 
    if_else(str_detect(icd, "401|402|403|404|405"), "HyperTn",
    if_else(str_detect(icd, paste("4", 10:14, sep = "", collapse = "|")), "IHD",
    if_else(str_detect(icd, paste("4", 15:17, sep = "", collapse = "|")), "Pul_HD",
    if_else(str_detect(icd, paste("4", 20:29, sep = "", collapse = "|")), "Other_HD",
    if_else(str_detect(icd, paste("4", 30:38, sep = "", collapse = "|")), "CerebroVD",
    if_else(str_detect(icd, paste("4", 40:49, sep = "", collapse = "|")), "Arteries",
    if_else(str_detect(icd, paste("4", 51:59, sep = "", collapse = "|")), "Veins",
    "NOI"))))))))
    
Potentially due to the data transferral from SPSS file to R, both ICD 9 and 10 codes are stored as factor instead of string in the variable HDIA3. To apply stringr for group, it is necessary to turn them into nice looking strings. 
First, there are many unused levels, which are categories that are actually not occupied by any individuals. The droplevels.factor commend can remove the categories that are without individuals. Next, the variable is treated as string. However, a little space is left at the end of all ICD codes. The str_trim is therefore used to remove the space at the end of all string ( str_trim(x,"right") ). 

Now we can proceed to group the ICD codes. The key issue here is to identify the pattern. ICD 9 are purely numbers, while ICD 10 starts with alphabets. str_detect is to find pattern that match and return the logical variable - TRUE or FALSE. So, to interpret the first line, using str_detect to find if there is any string in the variable icd that match 401 OR 402 OR ... . If there are indeed some of them match (return as TRUE by str_detect), please turn them as "HyperTn" in the new variable group_h. 
The use of sign | is to tell str_detect that there are few patterns to be matched and if any one of the given patterns match, please return TRUE. 
Knowing this format str_detect(variable with string, "pattern of string"), we can then pass patterns needed to group ICD. The second line is then use paste function to help create patterns that are separated by |. So, it means stick a 4 with 10, 11, 12, 13, and 14. (10:14 means from 10 to 14). The separation between 4 and 10, 4 and 11 ...  is nothing. All numbers are put together in one string and cut by the sign |. 

With stringr and the base function, the ICD codes are quickly organised into groups needed for later modelling. 
  
  