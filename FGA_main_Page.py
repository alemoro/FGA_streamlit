import streamlit as st
from PIL import Image

icon = Image.open("./img/Neurospector.png")
st.set_page_config(
    page_title="FGA data processing and analysis",
    page_icon=icon,
)

with st.sidebar:
    logo = Image.open("./img/FGA_Neurospector.png")
    st.image(logo)

st.write("# Welcome to the FGA data processing and analysis page!")

st.markdown("""
            This is were I should write some more info.

            ### What can you do?
            - Multidimensional analysis
            - ELISA analysis
            - Interact with SQLite database from CellProfiler

            ### Any other ideas?
            I don't have the time to do them!
            """)
