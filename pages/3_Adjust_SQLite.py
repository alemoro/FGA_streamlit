import streamlit as st
from PIL import Image

icon = Image.open("./img/Neurospector.png")

st.set_page_config(
    page_title="FGA ELISA analysis",
    page_icon=icon,
    layout='wide',
)
 

# Add a quick start on the sidebar
with st.sidebar:
    logo = Image.open("./img/FGA_Neurospector.png")
    st.image(logo)

st.subheader('How to do it', divider=True)
st.markdown("""
            OK, so here is the thing. I am about to leave the FGA, and I don't know how to make this app work with SQLite databases. I tried for a few hours, and I don't see how, so I will teach you how to do it in python.  
            My recomendation is to install python (if you don't already have it) via *anaconda*.  
            Once you installed it, you can open a Jupyter workbook (a bit easier to work with), and follow my instructions.  
            """)              
st.markdown("""
            **Load modules**  
            The first thing that we need to do is to load the modules that we will use:
            """)
st.code("""
        import sqlite3
        from tkinter import Tk
        from tkinter.filedialog import askopenfilename
        """)
st.markdown("""
            They should all be installed by default with anaconda, if not you should be able to find the explanation online.  
            The next step is to locate the SQLite database from CellProfiler and open it
            """)
st.code("""
        # Open file selection dialog to choose the SQLite database file
        Tk().withdraw()
        db_file = askopenfilename(title="Select SQLite Database File", filetypes=[("SQLite databases", "*.db")])
        # Open SQLite database
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        """)
st.markdown("""
            One of the things that you would like to do is to add a new column with metadata information, such as *genotype*, *treatment* or similar.  
            To do so we need to first create a new column in the database, then to find add a link between an existing column (*i.e.* the **WellID**) like in the examples below.
            """)
st.code('''
        # Loop through the tables that contain "_Per_Image" in the name
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = cursor.fetchall()
        for table in tables:
            if "Per_Image" in table[0]:
                # Add the new columns to the table
                cursor.execute(f"ALTER TABLE {table[0]} ADD COLUMN Genotype")
                cursor.execute(f"ALTER TABLE {table[0]} ADD COLUMN Treatment")
        
            # Update the "Genotype" column based on the condition
            cursor.execute(f"""
                UPDATE {table[0]}
                SET Genotype = 
                    CASE
                        WHEN Image_Metadata_WellID LIKE '%3' THEN 'Control'
                        WHEN Image_Metadata_WellID LIKE '%6' THEN 'Control'
                        WHEN Image_Metadata_WellID LIKE '%9' THEN 'Control'
                        WHEN Image_Metadata_WellID LIKE '%4' THEN 'D207G'
                        WHEN Image_Metadata_WellID LIKE '%7' THEN 'D207G'
                        WHEN Image_Metadata_WellID LIKE '%10' THEN 'D207G'
                        ELSE 'S241fs'
                    END
            """)
        
            # Update the "Treatment" column based on the condition
            cursor.execute(f"""
                UPDATE {table[0]}
                SET Treatment = 
                    CASE
                        WHEN Image_Metadata_WellID LIKE 'A%' THEN 'Untreated'
                        WHEN Image_Metadata_WellID LIKE 'B%' THEN 'Untreated'
                        WHEN Image_Metadata_WellID LIKE 'C%' THEN 'Lenti'
                        WHEN Image_Metadata_WellID LIKE 'D%' THEN 'Lenti'
                        ELSE 'ASO'
                    END
            """)

            # Save the table changes to the database
            conn.commit()
        ''')
st.markdown("""
            In the example we created two new variables:
            - Genotype: in this case we were assigning the name of the genotype based on the **column** position on the plate (*%3*)
            - Treatment: in this case we were assigning the name of the genotype based on the **row** position on the plate (*A%*)
              
            Additionally, you might want to chage the location of your images. Since CellProfiler add a link to the images from when it was first run, and you might have the need to move the files somewhere else.  
            This is simply done with those lines of code:
            """)
st.code('''
        # Loop through the tables that contain "_Per_Image" in the name
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = cursor.fetchall()
        for table in tables:
            if "Per_Image" in table[0]:
                cursor.execute(f"PRAGMA table_info({table[0]})")
                columns = cursor.fetchall()
                image_columns = [column[1] for column in columns if any(name in column[1] for name in "Image_PathName_")]
                for column in image_columns:
                    cursor.execute(f"UPDATE {table[0]} SET {column} = REPLACE({column}, '/scistor/', '//scistor.vu.nl/shares/')")
                    cursor.execute(f"UPDATE {table[0]} SET {column} = REPLACE({column}, '/', '\\')")    
                cursor.connection.commit()
        ''')
st.markdown("""
            In the example we changed the patch in two steps:
            - 1: going from "*/scistor/*" to "*//scistor.vu.nl/shares/*". This is important if the data was analyzed via BAZIS, since it use a different path to SciStor compared to Windows/MAC.
            - 2: change the file separators from "/" to "\\". This seems trivial and not important, but it is actually essential, since otherwise Windows will not recognize the name a valid file path. If you use a MAC it should not be needed, but you should double check.
              
            If you modify the location of you data you also have to modify the path in the *.properties file that CellProfiler created. Which should look something like this:
            """)
st.code('''
        #Tue May  2 08:49:16 2023
        # ==============================================
        #
        # CellProfiler Analyst 3.0 properties file
        #
        # ==============================================

        # ==== Database Info ====
        db_type         = sqlite
        db_sqlite_file  = \\scistor.vu.nl\shares\BETA-NeuroSciences-FGA-Screen\External\YourFolder\Your_CellProfiler_Analysis_Database.bd
        ''')
