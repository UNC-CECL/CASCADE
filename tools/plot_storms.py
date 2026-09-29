import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio

# get the plotly plots to show up
pio.renderers.default = 'browser'

# load the data
storms = pd.read_csv(r"C:\Users\agfig\model\calibration\storms_berm_2pt0_slope_0pt06.csv")
years = storms.time.values
durs = storms.duration.values
storms_per_yr, bins = np.histogram(years, bins=20)

# average storms
avg_storms = np.mean(storms_per_yr)
std_storms = np.std(storms_per_yr)
print(avg_storms, std_storms)

# default histogram
bins = int(max(durs))+1
fig = px.histogram(storms, x="duration", nbins=bins)

# updates the figure parameters to have an outline
fig.update_traces(marker_color='cornflowerblue', marker_line_color='black',
                  marker_line_width=1.0, opacity=1)
fig.update_layout(
    font=dict(
        size=18
    ),
    xaxis = dict(
        tickmode = 'linear',
        dtick = 2
    ),
    xaxis_tickangle=-90,
    xaxis_title="duration (hours)",
    yaxis_title="number of storms"
)
fig.show()

# default histogram
bins = int(max(durs))+1
fig = px.histogram(storms, x="time", nbins=bins)
vals = np.arange(2004, 2025, 1)
labels = []
for v in vals:
    labels.append(str(v))

# updates the figure parameters to have an outline
fig.update_traces(marker_color='cornflowerblue', marker_line_color='black',
                  marker_line_width=1.0, opacity=1)
fig.update_layout(
    font=dict(
        size=18
    ),
    xaxis = dict(
        tickmode = 'array',
        tickvals = np.arange(1,22,1),
        ticktext = labels
    ),
    xaxis_tickangle=-90,
    xaxis_title="year",
    yaxis_title="number of high water events"
)
fig.show()


# plot duration of each storm
# plot color by "time"
fig = px.scatter(
    storms,
    x=storms.index,
    y="duration",
    color=storms.time,
    color_continuous_scale = px.colors.cyclical.mygbm,
    # color=storms.time.astype(str),  # plotly automatically uses continuous colors bars
    # # for numbers and discrete color bars for strings
    # color_discrete_sequence=px.colors.qualitative.Plotly,
    symbol="time",
)
fig.update_traces(marker_size=10)
fig.update_layout(
    font=dict(
        size=14
    ),
    xaxis_title="storm number",
    yaxis_title="duration (hours)",
    legend_title_text="model year",
    showlegend=False,
    # legend_orientation="h",
    coloraxis_colorbar=dict(title='model year')
)
# fig.update_xaxes(showgrid=False)
# fig.update_yaxes(showgrid=False)
fig.show()