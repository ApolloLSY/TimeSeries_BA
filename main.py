import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from biolearn.data_library import DataLibrary
from biolearn.model_gallery import ModelGallery
import statsmodels.api as sm
from scipy.stats import boxcox
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox

def save_acf_pacf_to_pdf(df, max_lags=20, filename='ACF_PACF_Plots.pdf'):
    with PdfPages(filename) as pdf:
        for index, series in df.iterrows():
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            plot_acf(series, lags=max_lags, ax=ax1)
            ax1.set_title(f'Autocorrelation Function for Site {index}')
            plot_pacf(series, lags=max_lags, ax=ax2)
            ax2.set_title(f'Partial Autocorrelation Function for Site {index}')
            pdf.savefig(fig)
            plt.close(fig)  

def plot_first_five_sites(df):
    fig, axes = plt.subplots(5, 2, figsize=(12, 20), dpi=100) 
    for i, (index, series) in enumerate(df.iloc[:5].iterrows()):
        plot_acf(series, lags=20, ax=axes[i, 0], title=f'Site {index} ACF')
        plot_pacf(series, lags=20, ax=axes[i, 1], title=f'Site {index} PACF')
    plt.tight_layout()
    plt.savefig('First_Five_Sites_ACF_PACF.png')
    plt.close()
    
def plot_observations_and_forecasts(forecast_df, stationary_df):
    plt.figure(figsize=(12, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, 5)) 
    for (index, row), color in zip(forecast_df.iloc[:5].iterrows(), colors):
        observations = stationary_df.loc[row['Index']]
        observation_times = range(1, len(observations) + 1)
        forecast_values = [row[f'Forecast_{90 + i}'] for i in range(1, 6)]
        forecast_times = range(len(observations) + 1, len(observations) + 6)
        plt.plot(observation_times, observations, color=color, label=f'Site {row["Index"]} Observations', linestyle='-')
        plt.plot(forecast_times, forecast_values, color=color, label=f'Site {row["Index"]} Forecast', linestyle='--')
    plt.title('Observations and Forecasts for the First Five Sites')
    plt.xlabel('Time Point')
    plt.ylabel('Value')
    plt.ylim(0, 1)  
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.grid(True)
    plt.savefig('Observations_Forecasts_First_Five_Sites.png')
    plt.show()
   
gallery = ModelGallery()
datasets = ["GSE41169", "GSE40279", "GSE51057", "GSE42861", "GSE51032", 
            "GSE73103", "GSE69270", "GSE64495", "GSE30870", "GSE52588"]

predictions = []

for i, dataset_name in enumerate(datasets):
    data_source = DataLibrary().get(dataset_name)
    exec(f"data{i+1} = data_source.load()")
    exec(f"horvath_results{i+1} = gallery.get('Horvathv2').predict(eval(f'data{i+1}'))")
    predictions.append(eval(f"horvath_results{i+1}"))

prediction=pd.concat(predictions)

prediction['Predicted'] = prediction['Predicted'].round().astype(int)

prediction.to_csv("E:/111/age_prediction.csv",index=True,header=True,sep="\t")

panel = pd.DataFrame(0, index=list(data1.dnam.index), columns=range(1, 101))
count = pd.DataFrame(0, index=['count'], columns=range(1, 101))

for i, pred in enumerate(prediction['Predicted']):
    if pred in panel.columns:
        index_name = prediction.index[i]
        for j, dataset_name in enumerate(datasets):
            data = globals()[f'data{j+1}']  
            if index_name in data.dnam.columns:
                panel[pred] = panel[pred].add(data.dnam[index_name])
                count[pred] += 1
    else:
        print("100 is not enough")

count_transposed = count.T

plt.figure(figsize=(10, 6))
plt.scatter(count_transposed.index, count_transposed['count'])
plt.title('Count of Predictions')
plt.xlabel('Prediction Index')
plt.ylabel('Count')
plt.grid(True)
plt.show()

zero_count_columns = count.columns[count.iloc[0] == 0]
print("Count columns with zero count:", zero_count_columns)

panel_weighted = panel.copy()
for col in panel.columns:
    if col in zero_count_columns:
        panel_weighted[col] = 0
    else:
        panel_weighted[col] = panel[col] / count[col].iloc[0]

cg_strings = [
     "cg12140144", "cg26933021", "cg20822990", "cg07312601", "cg09993145", "cg23605843", 
    "cg25410668", "cg17879376", "cg14962509", "cg24375409", "cg22851420", "cg24107728", 
    "cg14614643", "cg00257455", "cg23045908", "cg15201877", "cg18933331", "cg05675373", 
    "cg19269039", "cg16008966", "cg14565725", "cg05940231", "cg03984502", "cg25256723", 
    "cg16054275", "cg01459453", "cg16599143", "cg02275294", "cg21870884", "cg10501210", 
    "cg02901139", "cg11298786", "cg09809672", "cg05940691", "cg12869659", "cg23028740", 
    "cg10959651", "cg00522231", "cg01752203", "cg23643435", "cg01243072", "cg07589899", 
    "cg10855531", "cg22943590", "cg23398076", "cg12082609", "cg01447660", "cg02085953", 
    "cg15149655", "cg22809047", "cg06639320", "cg22454769", "cg00017842", "cg23606718", 
    "cg22061831", "cg12757011", "cg01620164", "cg12105450", "cg00760938", "cg10376763", 
    "cg23077820", "cg08166272", "cg10523019", "cg23462687", "cg20669012", "cg03183882", 
    "cg15910502", "cg12941369", "cg02244028", "cg24888989", "cg00702638", "cg07303143", 
    "cg26614073", "cg15988232", "cg16933388", "cg19381811", "cg03019000", "cg01844642", 
    "cg04474832", "cg03891319", "cg03607117", "cg22264409", "cg11205552", "cg17321954", 
    "cg06796779", "cg18303397", "cg09025210", "cg14423778", "cg15277914", "cg07553761", 
    "cg06737494", "cg01059398", "cg26824216", "cg25478614", "cg07110949", "cg23239150", 
    "cg21254939", "cg23995914", "cg23836737", "cg05960024", "cg10699857", "cg05024939", 
    "cg05106770", "cg06690548", "cg02650266", "cg01511232", "cg25148589", "cg24843443", 
    "cg03364683", "cg21815258", "cg12608692", "cg20755989", "cg12238343", "cg02328239", 
    "cg21878650", "cg26921969", "cg10837404", "cg01883408", "cg06448705", "cg16983159", 
    "cg08234504", "cg11006267", "cg21874213", "cg23500537", "cg26843711", "cg08587542", 
    "cg10345936", "cg16281600", "cg03555227", "cg14345676", "cg14314729", "cg23517605", 
    "cg01570885", "cg23375552", "cg20052760", "cg16867657", "cg21572722", "cg00194146", 
    "cg01527307", "cg22736354", "cg10699171", "cg06493994", "cg04424621", "cg02281167", 
    "cg03771840", "cg06685111", "cg06462220", "cg08420066", "cg21467614", "cg12753631", 
    "cg18501647", "cg04576021", "cg10192196", "cg18468088", "cg03894990", "cg01740766", 
    "cg16255583", "cg04642300", "cg17266282", "cg07095347", "cg00073460", "cg16333846", 
    "cg13221458", "cg05468948", "cg00795927", "cg08911208", "cg22372849", "cg16012294", 
    "cg22679120", "cg13931228", "cg27009703", "cg11671968", "cg26312920", "cg19663246",
    "cg14396995", "cg18442362", "cg09748749", "cg20692569", "cg23857078",
    "cg21743182", "cg09436502", "cg00503840", "cg04084157", "cg14175438", "cg20665157",
    "cg21184711", "cg02383785", "cg04528819", "cg20426994", "cg08097417", "cg02821342",
    "cg20397034", "cg03473532", "cg08280936", "cg08540945", "cg18769120", "cg26101086",
    "cg19859445", "cg07502389", "cg18267374", "cg00582628", "cg16419235", "cg23710218",
    "cg07583137", "cg12402251", "cg19497517", "cg13586038", "cg19724470", "cg07211259",
    "cg07158339", "cg24046474", "cg14059835", "cg10570177", "cg13649056", "cg13734401",
    "cg26581729", "cg06231995", "cg14411282", "cg12530994", "cg23754392", "cg06908778",
    "cg22796704", "cg01560871", "cg04268405", "cg18738190", "cg04126866", "cg25427880",
    "cg09671951", "cg06888746", "cg24838825", "cg13848598", "cg07906193", "cg12776156",
    "cg05928581", "cg17627559", "cg23091758", "cg04940570", "cg10825530", "cg20654468",
    "cg26552743", "cg21992250", "cg15015340", "cg22843803", "cg02532488", "cg13547237",
    "cg06419846", "cg05496363", "cg20063906", "cg12328429", "cg25969122", "cg18633600",
    "cg08622677", "cg01820374", "cg25719851", "cg13828440", "cg26986871", "cg00748589",
    "cg21747310", "cg00431549", "cg13302154", "cg13909661", "cg19722847", "cg26311454",
    "cg08900043", "cg00753885", "cg18573383", "cg15405572", "cg01528542", "cg08993878",
    "cg22827210", "cg03670162", "cg04596060", "cg21907579", "cg10281002", "cg07172885",
    "cg20404336", "cg18582260", "cg22179082", "cg06648759", "cg23357533", "cg20102280",
    "cg13767001", "cg23389651", "cg09646392", "cg00593462", "cg07115626", "cg06738602",
    "cg03032497", "cg18771300", "cg24058132", "cg13027206", "cg15480367", "cg12177001",
    "cg23709172", "cg14334310", "cg04875128", "cg21296230", "cg02071305", "cg03361973",
    "cg16717122", "cg04858164", "cg25005357", "cg08454546", "cg21801378", "cg15188939",
    "cg20540209", "cg05542681", "cg08949164", "cg09183146", "cg00454305", "cg02871659",
    "cg00991848", "cg02331561", "cg26974444", "cg06112560", "cg10917602", "cg09155044",
    "cg04031656", "cg03746976", "cg18693704", "cg00658652", "cg03991512", "cg22947000",
    "cg07082267", "cg03486383", "cg02228185", "cg23668631", "cg14522800", "cg25135555",
    "cg13029847", "cg06144905", "cg11896923", "cg06874016", "cg25809905", "cg04267101",
    "cg22507023", "cg02867102", "cg13093111", "cg13683374", "cg10644544", "cg03643998",
    "cg11620135", "cg25436157", "cg24217948", "cg17589341", "cg17243289", "cg12459502",
    "cg19283806", "cg10052840", "cg26005082", "cg11766468", "cg10586358", "cg14556683",
    "cg26842024", "cg18335931", "cg17861230", "cg10498798", "cg27212234", "cg17110586",
    "cg21944491", "cg21940708", "cg15811427", "cg06458239", "cg10729426", "cg21911021",
    "cg24834740", "cg19702785", "cg07547549", "cg12303084", "cg27544190", "cg01262913",
    "cg17274064", "cg09428349", "cg10636297", "cg12373771", "cg05442902", "cg16612562",
    "cg19015086", "cg21205978", "cg01949403", "cg08415592", "cg19853760", "cg23124451",
    "cg25459323", "cg27187881", "cg00343092", "cg00347775", "cg27131176", "cg22982767",
    "cg07016730", "cg01892695"
]

core_cpg = pd.DataFrame(cg_strings, columns=['CpG'])

Pseudo_TimeSeries_100 = panel_weighted[panel_weighted.index.isin(core_cpg['CpG'])]
Pseudo_TimeSeries_90=Pseudo_TimeSeries_100.iloc[:,:90]
stats = Pseudo_TimeSeries_90.describe()

stats_T=stats.T
stats_T.to_csv("E:/111/Stats_Transformed.csv",index=True,header=True,sep="\t")
Pseudo_TimeSeries_90.to_csv("E:/111/Pseudo_TimeSeries_90.csv",index=True,header=True,sep="\t")

df = Pseudo_TimeSeries_90
total_columns = df.shape[1]

nan_count_per_row = df.isnull().sum(axis=1)
nan_percentage_per_row = nan_count_per_row / total_columns

filtered_df = df[nan_percentage_per_row <= 0.05]
df_filtered_filled = filtered_df.fillna(method='ffill', axis=1)
df_filtered_filled.iloc[:, 9] = np.nan
df_filtered_filled.iloc[:, 9]=df.interpolate(method='linear', axis=1)

mean_values = df_filtered_filled.mean(axis=1)  
std_values = df_filtered_filled.std(axis=1)    
cutoff = std_values * 4
lower = mean_values - cutoff
upper = mean_values + cutoff

outliers = (df_filtered_filled.lt(lower, axis=0) | df_filtered_filled.gt(upper, axis=0))

df_exclude_outliers =df_filtered_filled.copy()
df_exclude_outliers[outliers] = np.nan
df_exclude_outliers=df_exclude_outliers.fillna(method='bfill', axis=1)

plt.figure(figsize=(18, 24))
num_rows = len(df_filtered_filled)
colors = plt.cm.viridis(np.linspace(0, 1, num_rows))
for (index, row), color in zip(df_filtered_filled.iterrows(), colors):
    plt.plot(df_filtered_filled.columns, row, label=index, color=color)
plt.title('Time Series of Each Site')
plt.xlabel('Year')
plt.ylabel('Methylation Value')
plt.show()

df_exclude_outliers.to_csv("E:/111/Pseudo_TimeSeries_no_outliers.csv",index=True,header=True,sep="\t")
df_ex_outliiers_demo=df_exclude_outliers.iloc[:10,:]

plt.figure(figsize=(18, 6))
for index, row in df_ex_outliiers_demo.iterrows():
    plt.plot(df_ex_outliiers_demo.columns, row, label=index)
plt.title('Time Series of Each Site')
plt.xlabel('Year')
plt.ylabel('Value')
plt.legend()
plt.show()

adf_results_1 = []
for index, row in df_exclude_outliers.iterrows():
    result = adfuller(row.dropna(), autolag='AIC')  
    adf_results_1.append({
        'Row Index': index,
        'ADF Statistic': result[0],
        'p-value': result[1],
        'Critical Value (1%)': result[4]['1%'],
        'Critical Value (5%)': result[4]['5%'],
        'Critical Value (10%)': result[4]['10%'],
        'Stationarity': "Stationary" if result[1] < 0.05 else "Non-stationary"
    })
adf_results_df_1 = pd.DataFrame(adf_results_1)
adf_results_df_1.to_csv("E:/111/ADF_test_1.csv",index=True,header=True,sep="\t")

stationary_df_1 = df_exclude_outliers.loc[adf_results_df_1[adf_results_df_1['Stationarity'] == 'Stationary']['Row Index']]
non_stationary_df_1 = df_exclude_outliers.loc[adf_results_df_1[adf_results_df_1['Stationarity'] == 'Non-stationary']['Row Index']]

stationary_df_1.to_csv("E:/111/ADF_stationary_df_1.csv",index=True,header=True,sep="\t")
non_stationary_df_1.to_csv("E:/111/ADF_nonstationary_df_1.csv",index=True,header=True,sep="\t")

plt.figure(figsize=(18, 18))  
num_rows = len(stationary_df_1)
colors = plt.cm.viridis(np.linspace(0, 1, num_rows)) 
for (index, row), color in zip(stationary_df_1.iterrows(), colors):
    plt.plot(stationary_df_1.columns, row, label=str(index), color=color)  
plt.title('Time Series of Each Site')
plt.xlabel('Year')
plt.ylabel('Methylation Value')
plt.legend(title='Index', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

plt.figure(figsize=(18, 18))  
num_rows = len(non_stationary_df_1)
colors = plt.cm.viridis(np.linspace(0, 1, num_rows))  
for (index, row), color in zip(non_stationary_df_1.iterrows(), colors):
    plt.plot(non_stationary_df_1.columns, row, label=str(index), color=color)  
plt.title('Time Series of Each Site')
plt.xlabel('Year')
plt.ylabel('Methylation Value')
plt.show()

save_acf_pacf_to_pdf(stationary_df_1)
print("All ACF and PACF plots have been saved to ACF_PACF_Plots.pdf.")

plot_first_five_sites(stationary_df_1)
print("ACF and PACF plots for the first five sites have been saved to 'First_Five_Sites_ACF_PACF.png'.")

orders = [(1, 0), (1, 0), (1, 1), (2, 0), (1, 0),
          (1,0),(1,0),(1,0),(1,0),(1,0),
          (1,0),(1,1),(1,0),(1,0),(1,0),
          (1,0),(1,1),(1,0),(1,0),(1,0),
          (1,0),(1,1),(1,0),(1,1),(2,0),
          (1,1),(1,1),(2,0),(2,0),(2,0),
          (2,0),(1,0),(1,0),(1,1),(1,0),
          (1,0),(0,0),(1,0),(0,0),(1,0),
          (1,0),(1,0),(1,0),(1,0),(1,0),
          (1,0),(1,0),(1,0),(1,0),(1,0),
          (0,0),(1,0),(1,1),(1,0),(1,1),
          (1,0),(1,0),(1,1),(1,0),(1,0),
          (1,2),(1,0),(1,0)
          ]

results = []

for index, (p, q) in zip(stationary_df_1.index, orders):
    series = stationary_df_1.loc[index]
    model = ARIMA(series, order=(p, 0, q))
    result = model.fit()
    params = {
        'Index': index,
        'p': p,
        'q': q,
        'const': result.params.get('const', 0), 
        'sigma2': result.params.get('sigma2', 0) 
    }
    for i in range(1, p + 1):
        params[f'ar.L{i}'] = result.params.get(f'ar.L{i}', 0)  
    for i in range(1, q + 1):
        params[f'ma.L{i}'] = result.params.get(f'ma.L{i}', 0)  
    results.append(params)

results_df_1 = pd.DataFrame(results)
results_df_1.to_csv("E:/111/ARMA_65_result.csv",index=True,header=True,sep="\t")

error_results = []
pdf_pages = PdfPages('Time_Series_Diagnostics.pdf')

for index, (p, q) in zip(stationary_df_1.index, orders):
    series = stationary_df_1.loc[index]
    model = ARIMA(series, order=(p, 0, q))
    result = model.fit()
    residuals = result.resid
    fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    ax[0].plot(residuals)
    ax[0].set_title(f'Residuals for Series at Index {index}')
    sm.graphics.tsa.plot_acf(residuals, lags=20, ax=ax[1], title=f'ACF for Series at Index {index}')
    pdf_pages.savefig(fig)
    plt.close(fig)  
    lb_test = acorr_ljungbox(residuals, lags=[20], return_df=True)
    lb_pvalue = lb_test['lb_pvalue'].values[0]
    jb_test = sm.stats.jarque_bera(residuals)
    jb_pvalue = jb_test[1]
    error_results.append({
        'Index': index,
        'p': p,
        'q': q,
        'LB p-value': lb_pvalue,
        'Mean of Residuals': residuals.mean(),
        'Std of Residuals': residuals.std()
    })
pdf_pages.close()

error_results_df = pd.DataFrame(error_results)
error_results_df.to_csv("E:/111/63_error_results_df.csv",index=True,header=True,sep="\t")

forecast_results = []

for index, (p, q) in zip(stationary_df_1.index, orders):
    series = stationary_df_1.loc[index]
    model = ARIMA(series, order=(p, 0, q))
    fitted_model = model.fit()
    
    
    forecast = fitted_model.get_forecast(steps=5)
    mean_forecast = forecast.predicted_mean
    conf_int = forecast.conf_int()
    
    
    forecast_dict = {
        'Index': index,
    }
    
    for i, value in enumerate(mean_forecast, start=1):
        forecast_dict[f'Forecast_{90 + i}'] = value
        forecast_dict[f'CI_Lower_{90 + i}'] = conf_int.iloc[i - 1, 0]
        forecast_dict[f'CI_Upper_{90 + i}'] = conf_int.iloc[i - 1, 1]
    
    forecast_results.append(forecast_dict)

forecast_df = pd.DataFrame(forecast_results)
forecast_df.to_csv("E:/111/63_forecast.csv",index=True,header=True,sep="\t")

plot_observations_and_forecasts(forecast_df, stationary_df_1)

data_boxcox_diff1 = pd.DataFrame()

results_2 = [] 
for index, row in non_stationary_df_1.iterrows():
    data, lambda_ = boxcox(row[row > 0])
    data_diff = np.diff(data, n=1) 
    temp_df = pd.DataFrame(data_diff.reshape(1, -1))
    data_boxcox_diff1 = pd.concat([data_boxcox_diff1, temp_df], ignore_index=True)
    adf_result = adfuller(data_diff)
    results_2.append({
        "index": index,
        "lambda": lambda_,
        "ADF Statistic": adf_result[0],
        "p-value": adf_result[1],
        'Critical Value (1%)': adf_result[4]['1%'],
        'Critical Value (5%)': adf_result[4]['5%'],
        'Critical Value (10%)': adf_result[4]['10%'],
        'Stationarity': "Stationary" if adf_result[1] < 0.05 else "Non-stationary"
    })
results_df_2 = pd.DataFrame(results)
data_boxcox_diff1.index = non_stationary_df_1.index

data_boxcox_diff1.to_csv("E:/111/data_boxcox_diff1.csv",index=True,header=True,sep="\t")
results_df_2.to_csv("E:/111/ADF_test_2.csv",index=True,header=True,sep="\t")

stationary_df_2 =  data_boxcox_diff1.loc[results_df_2[results_df_2['Stationarity'] == 'Stationary']['index']]
non_stationary_df_2 = data_boxcox_diff1.loc[results_df_2[results_df_2['Stationarity'] == 'Non-stationary']['index']]

plt.figure(figsize=(18, 18))  
num_rows = len(stationary_df_2)
colors = plt.cm.viridis(np.linspace(0, 1, num_rows))  
for (index, row), color in zip(stationary_df_2.iterrows(), colors):
    plt.plot(stationary_df_2.columns, row, label=str(index), color=color)  
plt.title('Time Series of Each Site')
plt.xlabel('Year')
plt.ylabel('Methylation Value')
plt.show()

plt.figure(figsize=(18, 12)) 
num_rows = len(non_stationary_df_2)
colors = plt.cm.viridis(np.linspace(0, 1, num_rows))  
for (index, row), color in zip(non_stationary_df_2.iterrows(), colors):
    plt.plot(non_stationary_df_2.columns, row, label=str(index), color=color)  
plt.title('Time Series of Each Site')
plt.xlabel('Year')
plt.ylabel('Methylation Value')
plt.show()

save_acf_pacf_to_pdf(stationary_df_2)
print("All ACF and PACF plots have been saved to ACF_PACF_Plots.pdf.")
plot_first_five_sites(stationary_df_2)
print("ACF and PACF plots for the first five sites have been saved to 'First_Five_Sites_ACF_PACF.png'.")

orders = [(1, 1), (0, 0), (0, 0), (1, 1), (0, 0),
          (1,1),(1,1),(0,0),(2,3),(1,1),
          (1,1),(3,1),(1,1),(1,1),(1,1),
          (0,0),(0,0),(1,1),(1,1),(1,1),
          (0,0),(0, 0),(0, 0),(0, 0),(1,1),
          (2,1),(1,1),(1,1),(2,1),(2,2),
          (1,1),(1,1),(1,1),(1,1),(2,1),
          (1,1),(1,1),(1,1),(1,1),(1,1),
          (0,0),(1,1),(2,1),(1,1),(0,0),
          (0,0),(0,0),(2,1),(3,2),(1,1),
          (1,1),(1,1),(1,1),(1,1),(1,1),
          (1,1),(0,0),(0,0),(1,1),(1,1),
          (0,0),(1,1),(1,1),(0,0),(1,1),
          (2,1),(0,0),(2,1),(0,0),(2,2),
          (1,1),(0,0),(1,1),(1,1),(1,1),
          (1,1),(0,0),(1,1),(1,1),(1,1),
          (0,0),(2,1),(0,0),(1,1),(0,0),
          (0,0),(0,0),(0,0),(1,1),(0,0),
          (2,1),(2,1),(0,0),(0,0),(0,0),
          (0,0),(0,0),(1,1),(2,1),(0,0),
          (1,1),(0,0),(1,1),(0,0),(1,1),
          (1,1),(1,1),(2,1),(1,1),(1,1),
          (1,1),(0,0),(0,0),(2,1),(0,0),
          (1,1),(2,1),(1,1),(2,1),(0,0),
          (2,1),(0,0),(0,0),(1,1),(2,1),
          (0,0),(0,0),(0,0),(1,1),(1,1),
          (1,1),(1,1),(1,1),(0,0),(1,1),
          (0,0),(1,1),(1,1),(0,0),(1,1),
          (0,0),(1,1),(1,1),(1,1),(1,1),
          (1,1),(2,1),(1,1),(2,1),(0,0),
          (1,1),(0,0),(1,1),(0,0),(0,0),
          (0,0),(2,1),(2,1),(1,1),(1,1),
          (0,0),(2,1),(2,1),(1,1),(2,3),
          (1,1),(1,1),(0,0),(1,1),(0,0),
          (0,0),(1,1),(1,1),(0,0),(1,1),
          (1,1),(0,0),(2,1),(1,1),(0,0),
          (0,0),(1,1),(0,0),(0,0),(1,1),
          (0,0),(1,1),(1,1),(1,1),(2,1),
          (1,1),(0,0),(1,1),(0,0),(1,1),
          (0,0),(1,1),(0,0),(1,1),(1,1),
          (0,0),(0,0),(0,0),(1,1),(0,0),
          (0,0),(1,1),(1,1),(2,1),(1,1),
          (0,0),(1,1),(1,1),(2,1),(2,1),
          (0,0),(0,0),(1,1),(1,1),(0,0),
          (2,2),(0,0),(2,1)
          ]

results = []
for index, (p, q) in zip(stationary_df_2.index, orders):
    series = stationary_df_2.loc[index]
    model = ARIMA(series, order=(p, 1, q))
    result = model.fit()
    params = {
        'Index': index,
        'p': p,
        'q': q,
        'const': result.params.get('const', 0),  
        'sigma2': result.params.get('sigma2', 0)  
    }
    for i in range(1, p + 1):
        params[f'ar.L{i}'] = result.params.get(f'ar.L{i}', 0) 
        
    for i in range(1, q + 1):
        params[f'ma.L{i}'] = result.params.get(f'ma.L{i}', 0) 
    results.append(params)

results_df_2 = pd.DataFrame(results)

forecast_results = []

for index, (p, q) in zip(stationary_df_2.index, orders):
    series = stationary_df_2.loc[index]
    model = ARIMA(series, order=(p, 0, q))
    fitted_model = model.fit()
    
    forecast = fitted_model.get_forecast(steps=5)
    mean_forecast = forecast.predicted_mean
    conf_int = forecast.conf_int()
    
    forecast_dict = {
        'Index': index,
    }

    for i, value in enumerate(mean_forecast, start=1):
        forecast_dict[f'Forecast_{90 + i}'] = value
        forecast_dict[f'CI_Lower_{90 + i}'] = conf_int.iloc[i - 1, 0]
        forecast_dict[f'CI_Upper_{90 + i}'] = conf_int.iloc[i - 1, 1]
    
    forecast_results.append(forecast_dict)
forecast_df = pd.DataFrame(forecast_results)
forecast_df.to_csv("E:/111/224_forecast.csv",index=True,header=True,sep="\t")
plt.show()

plot_observations_and_forecasts(forecast_df, stationary_df_2)
results_df_2.to_csv("E:/111/224_result.csv",index=True,header=True,sep="\t")





