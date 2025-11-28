import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

import doubleml as dml
from doubleml.datasets import make_confounded_plr_data, fetch_401K
import matplotlib.pyplot as plt
import plotly.io as pio
from sklearn.preprocessing import StandardScaler
import re


flux=pd.read_csv("data/key_sp_flux/key_sp_flux.csv",index_col=0)
selected=pd.read_csv("data/res_grow_wd_pFBA/selected_sample.csv")
flux=flux[flux.index.isin(selected.x)]
meta=pd.read_csv("data/test_cohort/test_16s/SraRunTable.csv",index_col=0)
meta["Library Name"]=meta["Library Name"].str.replace("[-.]","_",regex=True)
meta=meta[meta['Sample Name'].isin(selected.x)]
flux=flux.loc[meta["Library Name"].tolist()]
index=list()
for i in range(flux.shape[1]):
  tmp=np.mean(flux.iloc[:,i]!=0)
  if len(flux.iloc[:,i].unique())>3:
    index.append(i)
#flux = flux.drop(index="SZAXPI015258_57")
#meta=meta[meta["sample_id"]!="SZAXPI015258_57"]
flux=flux.iloc[:,index]
scaler = StandardScaler()


cf_y = 0.1
cf_d = 0.1
theta = float(0)

np.random.seed(42)
df = flux
meta.index=meta["Library Name"]
df["group"]=meta["diastolic_bp"]
df["age"]=meta["AGE"]
df["Adiponectin"]=meta["Adiponectin"]
df["BMI"]=meta["BMI"]
df["Body_Fat_Percentage"]=meta["Body_Fat_Percentage"]
tmp=meta.loc[:,['Glucose', 'Hb1Ac', 'HDL','hs-CRP','Insulin','LDL','Total_Cholesterol', 'Triglycerides', 'VLDL', 'Waist_circumference']]
df=pd.concat([df,tmp],axis=1)
scale_flux = scaler.fit_transform(df)
df=pd.DataFrame(scale_flux, columns=df.columns,index=df.index)
df=pd.concat([df,meta.loc[:,['city','sex', 'smoker','Medicament_use']]],axis=1)
df.iloc[:,20:30].hist(bins=30)
plt.show()
plt.close()
df.iloc[:, 1:30].plot(kind='box')
plt.show()
plt.close()
res=pd.DataFrame()
for i in range(len(df.columns)-19):
  d=str(df.columns[i])
  df1=df.loc[:,[d,'group', 'age', 'Adiponectin', 'BMI', 'Body_Fat_Percentage', 'Glucose', 'Hb1Ac',
  'HDL', 'hs-CRP', 'Insulin', 'Total_Cholesterol', 'Triglycerides', 'VLDL', 
  'Waist_circumference', 'city', 'sex', 'smoker', 'Medicament_use']]
  df1=pd.get_dummies(df1,columns=['city','Medicament_use', 'sex', 'smoker'])
  dml_data = dml.DoubleMLData(df1, 'group', d)
  print(dml_data)
  np.random.seed(123)
  dml_obj = dml.DoubleMLPLR(dml_data,
                          ml_l=RandomForestRegressor(),
                          ml_m=RandomForestRegressor(),
                          n_folds=5,
                          score='partialling out')
  dml_obj.fit()
  print(dml_obj)
  dml_obj.sensitivity_analysis(cf_y=cf_y, cf_d=cf_d, null_hypothesis=theta)
  print(dml_obj.sensitivity_summary)
  params=dml_obj.sensitivity_params
  tmp=pd.DataFrame({"ATE":dml_obj.coef,"CI2.5":dml_obj.confint().iloc[0,0],"CI97.5":dml_obj.confint().iloc[0,1],"pvalue":dml_obj.pval,
  "theta lower":params["theta"]["lower"],"theta upper":params["theta"]["upper"],"RV":params["rv"],"rva":params["rva"],
  "rmse ml_l":dml_obj.nuisance_loss["ml_l"][0],"rmse ml_m":dml_obj.nuisance_loss["ml_m"][0]})
  res=pd.concat([res,tmp],axis=0,ignore_index=True)
  fig=dml_obj.sensitivity_plot()
  pio.write_image(fig, "output/sensitivity/key_sp/"+d+"_sensitivity.png", scale=3)

res.index=df.columns[:-19]
res.to_csv("output/keysp_doubleml_diastolic.csv")


flux=pd.read_csv("data/key_sp_flux/key_sp_flux.csv",index_col=0)
selected=pd.read_csv("data/res_grow_wd_pFBA/selected_sample.csv")
flux=flux[flux.index.isin(selected.x)]
meta=pd.read_csv("data/test_cohort/test_16s/SraRunTable.csv",index_col=0)
meta["Library Name"]=meta["Library Name"].str.replace("[-.]","_",regex=True)
meta=meta[meta['Sample Name'].isin(selected.x)]
flux=flux.loc[meta["Library Name"].tolist()]
index=list()
for i in range(flux.shape[1]):
  tmp=np.mean(flux.iloc[:,i]!=0)
  if len(flux.iloc[:,i].unique())>3:
    index.append(i)
#flux = flux.drop(index="SZAXPI015258_57")
#meta=meta[meta["sample_id"]!="SZAXPI015258_57"]
flux=flux.iloc[:,index]
scaler = StandardScaler()


cf_y = 0.1
cf_d = 0.1
theta = float(0)

np.random.seed(42)
df = flux
meta.index=meta["Library Name"]
df["group"]=meta["Systolic_BP"]
df["age"]=meta["AGE"]
df["Adiponectin"]=meta["Adiponectin"]
df["BMI"]=meta["BMI"]
df["Body_Fat_Percentage"]=meta["Body_Fat_Percentage"]
tmp=meta.loc[:,['Glucose', 'Hb1Ac', 'HDL','hs-CRP','Insulin','LDL','Total_Cholesterol', 'Triglycerides', 'VLDL', 'Waist_circumference']]
df=pd.concat([df,tmp],axis=1)
scale_flux = scaler.fit_transform(df)
df=pd.DataFrame(scale_flux, columns=df.columns,index=df.index)
df=pd.concat([df,meta.loc[:,['city','sex', 'smoker','Medicament_use']]],axis=1)
df.iloc[:,20:30].hist(bins=30)
plt.show()
plt.close()
df.iloc[:, 1:30].plot(kind='box')
plt.show()
plt.close()
res=pd.DataFrame()
for i in range(len(df.columns)-19):
  d=str(df.columns[i])
  df1=df.loc[:,[d,'group', 'age', 'Adiponectin', 'BMI', 'Body_Fat_Percentage', 'Glucose', 'Hb1Ac',
  'HDL', 'hs-CRP', 'Insulin', 'Total_Cholesterol', 'Triglycerides', 'VLDL', 
  'Waist_circumference', 'city', 'sex', 'smoker', 'Medicament_use']]
  df1=pd.get_dummies(df1,columns=['city','Medicament_use', 'sex', 'smoker'])
  df1 = df1.dropna(subset=['group'])
  dml_data = dml.DoubleMLData(df1, 'group', d)
  print(dml_data)
  np.random.seed(123)
  dml_obj = dml.DoubleMLPLR(dml_data,
                          ml_l=RandomForestRegressor(),
                          ml_m=RandomForestRegressor(),
                          n_folds=5,
                          score='partialling out')
  dml_obj.fit()
  print(dml_obj)
  dml_obj.sensitivity_analysis(cf_y=cf_y, cf_d=cf_d, null_hypothesis=theta)
  print(dml_obj.sensitivity_summary)
  params=dml_obj.sensitivity_params
  tmp=pd.DataFrame({"ATE":dml_obj.coef,"CI2.5":dml_obj.confint().iloc[0,0],"CI97.5":dml_obj.confint().iloc[0,1],"pvalue":dml_obj.pval,
  "theta lower":params["theta"]["lower"],"theta upper":params["theta"]["upper"],"RV":params["rv"],"rva":params["rva"],
  "rmse ml_l":dml_obj.nuisance_loss["ml_l"][0],"rmse ml_m":dml_obj.nuisance_loss["ml_m"][0]})
  res=pd.concat([res,tmp],axis=0,ignore_index=True)
  fig=dml_obj.sensitivity_plot()
  pio.write_image(fig, "output/sensitivity/key_sp/"+d+"_sensitivity.png", scale=3)

res.index=df.columns[:-19]
res.to_csv("output/keysp_doubleml_Systolic_.csv")




### val ###


flux=pd.read_csv("data/key_sp_flux/val_key_sp_flux.csv",index_col=0)

meta=pd.read_csv("data/ZY_test/meta.csv",index_col=0)
meta=meta.loc[flux.index,:]
meta.BMI=pd.to_numeric(meta.BMI, errors='coerce')
meta=meta[meta.SBP.notna()]
meta=meta[meta.BMI.notna()]
flux=flux.loc[meta.index.tolist()]
index=list()
for i in range(flux.shape[1]):
  tmp=np.mean(flux.iloc[:,i]!=0)
  if len(flux.iloc[:,i].unique())>3:
    index.append(i)
#flux = flux.drop(index="SZAXPI015258_57")
#meta=meta[meta["sample_id"]!="SZAXPI015258_57"]
flux=flux.iloc[:,index]
scaler = StandardScaler()


cf_y = 0.1
cf_d = 0.1
theta = float(0)

np.random.seed(42)
df = flux
df["group"]=meta["DBP"]
df["age"]=meta.age
tmp=meta.loc[:,['BMI', 'HbA1C','Total.cholesterol', 'Triglycerides', 'HDL', 'LDL']]
df=pd.concat([df,tmp],axis=1)
scale_flux = scaler.fit_transform(df)
df=pd.DataFrame(scale_flux, columns=df.columns,index=df.index)
df=pd.concat([df,meta.loc[:,['sex', 'marriage']]],axis=1)
df.iloc[:,20:30].hist(bins=30)
plt.show()
plt.close()
df.iloc[:, 1:30].plot(kind='box')
plt.show()
plt.close()
res=pd.DataFrame()
for i in range(len(df.columns)-10):
  d=str(df.columns[i])
  df1=df.loc[:,[d,'group', 'age', 'BMI', 'HbA1C','Total.cholesterol', 'Triglycerides', 'HDL', 'LDL','sex','marriage']]
  df1=pd.get_dummies(df1,columns=['sex', 'marriage'])
  dml_data = dml.DoubleMLData(df1, 'group', d)
  print(dml_data)
  np.random.seed(123)
  dml_obj = dml.DoubleMLPLR(dml_data,
                          ml_l=RandomForestRegressor(),
                          ml_m=RandomForestRegressor(),
                          n_folds=5,
                          score='partialling out')
  dml_obj.fit()
  print(dml_obj)
  dml_obj.sensitivity_analysis(cf_y=cf_y, cf_d=cf_d, null_hypothesis=theta)
  print(dml_obj.sensitivity_summary)
  params=dml_obj.sensitivity_params
  tmp=pd.DataFrame({"ATE":dml_obj.coef,"CI2.5":dml_obj.confint().iloc[0,0],"CI97.5":dml_obj.confint().iloc[0,1],"pvalue":dml_obj.pval,
  "theta lower":params["theta"]["lower"],"theta upper":params["theta"]["upper"],"RV":params["rv"],"rva":params["rva"],
  "rmse ml_l":dml_obj.nuisance_loss["ml_l"][0],"rmse ml_m":dml_obj.nuisance_loss["ml_m"][0]})
  res=pd.concat([res,tmp],axis=0,ignore_index=True)
  fig=dml_obj.sensitivity_plot()
  pio.write_image(fig, "output/sensitivity/key_sp/"+d+"_sensitivity.png", scale=3)

res.index=df.columns[:-10]
res.to_csv("output/val_keysp_doubleml_diastolic.csv")





flux=pd.read_csv("data/key_sp_flux/val_key_sp_flux.csv",index_col=0)
meta=pd.read_csv("data/ZY_test/meta.csv",index_col=0)
meta=meta.loc[flux.index,:]
meta.BMI=pd.to_numeric(meta.BMI, errors='coerce')
meta=meta[meta.SBP.notna()]
meta=meta[meta.BMI.notna()]
flux=flux.loc[meta.index.tolist()]
index=list()
for i in range(flux.shape[1]):
  tmp=np.mean(flux.iloc[:,i]!=0)
  if len(flux.iloc[:,i].unique())>3:
    index.append(i)
#flux = flux.drop(index="SZAXPI015258_57")
#meta=meta[meta["sample_id"]!="SZAXPI015258_57"]
flux=flux.iloc[:,index]
scaler = StandardScaler()

cf_y = 0.1
cf_d = 0.1
theta = float(0)

np.random.seed(42)
df = flux
df["group"]=meta["SBP"]
df["age"]=meta.age
tmp=meta.loc[:,['BMI', 'HbA1C','Total.cholesterol', 'Triglycerides', 'HDL', 'LDL']]
df=pd.concat([df,tmp],axis=1)
scale_flux = scaler.fit_transform(df)
df=pd.DataFrame(scale_flux, columns=df.columns,index=df.index)
df=pd.concat([df,meta.loc[:,['sex', 'marriage']]],axis=1)
df.iloc[:,20:30].hist(bins=30)
plt.show()
plt.close()
df.iloc[:, 1:30].plot(kind='box')
plt.show()
plt.close()
res=pd.DataFrame()
for i in range(len(df.columns)-10):
  d=str(df.columns[i])
  df1=df.loc[:,[d,'group', 'age', 'BMI', 'HbA1C','Total.cholesterol', 'Triglycerides', 'HDL', 'LDL','sex','marriage']]
  df1=pd.get_dummies(df1,columns=['sex', 'marriage'])
  dml_data = dml.DoubleMLData(df1, 'group', d)
  print(dml_data)
  np.random.seed(123)
  dml_obj = dml.DoubleMLPLR(dml_data,
                          ml_l=RandomForestRegressor(),
                          ml_m=RandomForestRegressor(),
                          n_folds=5,
                          score='partialling out')
  dml_obj.fit()
  print(dml_obj)
  dml_obj.sensitivity_analysis(cf_y=cf_y, cf_d=cf_d, null_hypothesis=theta)
  print(dml_obj.sensitivity_summary)
  params=dml_obj.sensitivity_params
  tmp=pd.DataFrame({"ATE":dml_obj.coef,"CI2.5":dml_obj.confint().iloc[0,0],"CI97.5":dml_obj.confint().iloc[0,1],"pvalue":dml_obj.pval,
  "theta lower":params["theta"]["lower"],"theta upper":params["theta"]["upper"],"RV":params["rv"],"rva":params["rva"],
  "rmse ml_l":dml_obj.nuisance_loss["ml_l"][0],"rmse ml_m":dml_obj.nuisance_loss["ml_m"][0]})
  res=pd.concat([res,tmp],axis=0,ignore_index=True)
  fig=dml_obj.sensitivity_plot()
  pio.write_image(fig, "output/sensitivity/key_sp/"+d+"_sensitivity.png", scale=3)

res.index=df.columns[:-10]
res.to_csv("output/val_keysp_doubleml_systolic.csv")
