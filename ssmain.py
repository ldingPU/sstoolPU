import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt
from ssfunc import pf_calc
from ssfunc import lsm_sys_2droop
from ssfunc import lsm_sys_droop_gfl
from ssdata import case_3bus_2droop
from ssdata import case_3bus_droop_gfl
from itertools import product

st.set_page_config(layout="wide")
st.image('fig/TestSystem.png')
caption_html = """
<p style='
    text-align: center;
    color: black;          
    font-size: 18px;        
    font-family: Arial; 
    font-weight: bold;
'>
    System Configuration
</p>
"""
st.markdown(caption_html, unsafe_allow_html=True)


vector1 = ['GFM_Droop', 'GFL']
vector2 = ['GFM_Droop', 'GFL']
# Create all combinations and index them from 'com1' to 'com16'
combinations = list(product(vector1, vector2))
combination_named_index = {f'com{index + 1}': combination for index, combination in enumerate(combinations)}
selected_combination1 = combination_named_index['com1']
selected_combination2 = combination_named_index['com2']
# sidebar
sidebar1 = st.sidebar.selectbox(
    "What configuration do you want to select for the 1st generator?",
    ("GFM_Droop", "GFM_VSM", "GFL", "SG"), index=0, placeholder="Select configuration...",
)
sidebar2 = st.sidebar.selectbox(
    "What configuration do you want to select for the 2nd generator?",
    ("GFM_Droop", "GFM_VSM", "GFL", "SG"), index=0, placeholder="Select configuration...",
)
if (sidebar1,sidebar2) == ('GFM_Droop','GFM_Droop'):
    st.write('You selected GFM_Droop and GFM_Droop.')
    sysData = case_3bus_2droop()
    x, xdot = lsm_sys_2droop(sysData)
    Xss = [0.0000,0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,
           -0.1388,0.0386,-0.1930,0.5132,-0.1407,0.9993,0.5147,0.1411,
           0.1452,-0.0386,0.0844,-0.0002,0.5228,-0.1388,0.0386,-0.1930,
           0.5132,-0.1407,1.0249,-0.2814,0.5127,-0.1407,0.5127,-0.1407]
    stateVariableNames = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                          'theta2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                          'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']
elif (sidebar1,sidebar2) == ('GFM_Droop','GFL'):
    st.write('You selected GFM_Droop and GFL.')
    sysData = case_3bus_droop_gfl()
    x, xdot = lsm_sys_droop_gfl(sysData)
    Xss = [0.0000,0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,
           -0.1388,0.0386,-0.1930,0.5132,-0.1407,0.0000,0.0000,0.9993,
           0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,-0.1388,
           0.0386,-0.1930,0.5132,-0.1407,1.0249,-0.2814,0.5127,-0.1407,
           0.5127,-0.1407]
    stateVariableNames = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                          'theta2','epsilonL2','wf2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                          'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']
elif (sidebar1,sidebar2) == ('GFL','GFM_Droop'):
    st.write('You selected GFL and GFM_Droop.')
    sysData = case_3bus_droop_gfl()
    x, xdot = lsm_sys_droop_gfl(sysData)
    Xss = [0.0000,0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,
           -0.1388,0.0386,-0.1930,0.5132,-0.1407,0.0000,0.0000,0.9993,
           0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,-0.1388,
           0.0386,-0.1930,0.5132,-0.1407,1.0249,-0.2814,0.5127,-0.1407,
           0.5127,-0.1407]
    stateVariableNames = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                          'theta2','epsilonL2','wf2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                          'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']
else:
    st.write("The selected combination is not supported.")
Asys = xdot.jacobian(x)
for i in range(len(Xss)):
    Asys = Asys.subs([(x[i], Xss[i])])
Asys = np.array(Asys).astype(np.float64)
eigvals, eigenvectors = np.linalg.eig(Asys)
numeigs = len(eigvals)
lefteigenvectors = np.linalg.inv(eigenvectors)
pmatrix = np.multiply(eigenvectors,np.transpose(lefteigenvectors))
pmatrixabs = abs(pmatrix)

modeNames = ['mode{}'.format(i) for i in range(1,numeigs+1)]
# plot participation factor map
figheatmap = px.imshow(pmatrixabs,
                       labels=dict(x="Modes", y="State Variables"),
                       x = list(range(1,numeigs+1)),
                       y = stateVariableNames)
figheatmap.update_layout(height=800)

# Update font size, color and style
figheatmap.update_xaxes(title_font=dict(size=18,color='black',family='Arial'),tickfont=dict(size=14,color='black',family='Arial'))
figheatmap.update_yaxes(title_font=dict(size=18,color='black',family='Arial'),tickfont=dict(size=14,color='black',family='Arial'))

# Use text_input for manual number input
input_number = st.sidebar.text_input("Which mode do you want to select?", value='1') #(1-"+str(numeigs)+")
number = int(input_number)
# Check if the input is a number and within the desired range
if input_number:
    try:        
        number = int(input_number) # Convert input to an integer        
        if 1 <= number <= numeigs: # Check if the number is in the range
            df1 = pd.DataFrame(pmatrixabs, columns=modeNames)
            df1.insert(0, "statevariables", stateVariableNames, True)           
            df1.loc[df1[modeNames[number-1]] < 0.02, 'statevariables'] = 'Other states'  # Represent state variables with a relatively larger participation factor
            figpie1 = px.pie(df1, values=modeNames[number-1], names='statevariables', title='Participation factor analysis of mode '+str(number))
            figpie1.update_layout(title={'text':'Participation factor analysis of mode '+str(number),'x':0.415,'xanchor':'center'})
            df2 = pd.DataFrame(pmatrixabs, columns=modeNames)
            df2.insert(0, "statevariables", stateVariableNames, True)           
            df2.loc[df2[modeNames[number]] < 0.02, 'statevariables'] = 'Other states'  # Represent state variables with a relatively larger participation factor
            figpie2 = px.pie(df2, values=modeNames[number], names='statevariables', title='Participation factor analysis of mode '+str(number+1))
            figpie2.update_layout(title={'text':'Participation factor analysis of mode '+str(number+1),'x':0.415,'xanchor':'center'})

            # Update font size, color and style
            figpie1.update_layout(title={'text':'Participation factor analysis of mode '+str(number),'x':0.415,'xanchor':'center','font': {'size': 18, 'color': 'black','family': 'Arial'}}) # Update the layout for the title
            figpie1.update_traces(textfont={'size': 14, 'color': 'black','family': 'Arial'}) # Update the traces for the labels inside the pie chart
            figpie1.update_layout(legend_title_font={'size': 14, 'color': 'black'}, legend_font={'size': 14, 'color': 'black','family': 'Arial'}) # Update the legend font size and color      
            figpie2.update_layout(title={'text':'Participation factor analysis of mode '+str(number+1),'x':0.415,'xanchor':'center','font': {'size': 18, 'color': 'black','family': 'Arial'}}) # Update the layout for the title
            figpie2.update_layout(title={'text':'Participation factor analysis of mode '+str(number+1),'x':0.415,'xanchor':'center'})
            figpie2.update_traces(textfont={'size': 14, 'color': 'black','family': 'Arial'}) # Update the traces for the labels inside the pie chart
            figpie2.update_layout(legend_title_font={'size': 14, 'color': 'black'}, legend_font={'size': 14, 'color': 'black','family': 'Arial'}) # Update the legend font size and color            
        
            # Update the figure/chart position
            #figpie1.update_layout(margin=dict(t=10)) # Adjust left, right, top, bottom margins (l=xx, r=xx, t=xx, b=xx)
        
        else:
            st.error('Number out of range. Please enter a number between 1 and '+str(numeigs)+'.')
    except ValueError:        
        st.error('Invalid input. Please enter a number.') # Handle the case where input is not a number

col1, col2 = st.columns(2,gap="small")
with col1:
    st.plotly_chart(figheatmap, height=1200, theme="streamlit",use_container_width=True)
with col2:
    st.plotly_chart(figpie1, height=400, theme="streamlit",use_container_width=True)
    st.plotly_chart(figpie2, height=400, theme="streamlit",use_container_width=True)
    #eigvalsi = eigvals[number-1]  
    #st.text("real: "+str(eigvalsi.real))
    #st.text("imag: "+str(eigvalsi.imag))
    #st.text("freq: "+str(eigvalsi.imag/2/math.pi)+" Hz")
    #st.text("damp: "+str(-eigvalsi.real/np.sqrt(eigvalsi.real*eigvalsi.real+eigvalsi.imag*eigvalsi.imag)))

# plot table
colnew1, colnew2, colnew3 = st.columns(3,gap="small")
with colnew1:
    st.text("")
with colnew2:
    mode = range(1,len(eigvals)+1)
    realpart = eigvals.real
    imagpart = eigvals.imag
    frequency = eigvals.imag/2/math.pi
    dampingratio = -eigvals.real/np.sqrt(realpart*realpart+imagpart*imagpart)
    list_of_tuples = list(zip(mode, realpart, imagpart, frequency, dampingratio)) 
    df = pd.DataFrame(list_of_tuples,
                      columns = ["mode", "real", "imag", "freq(Hz)", "damping ratio"])
    st.table(df)
with colnew3:
    st.text("")

