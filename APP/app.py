# app.py

import streamlit as st
from Moving_Load_Opensees import (
    define_model_parameters,
    create_opensees_model,
    apply_moving_load,
    run_analysis,
    post_processing_results,
    plot_results
)

"""
Author: Mohammad Talebi-Kalaleh (talebika@ualberta.ca)
Feel free to reach out with any inquiries or suggestions.
"""

def main():
    # -- Title and Intro
    st.title("Moving Load Analysis in OpenSees")
    st.markdown(
        """
        **Author**: [Mohammad Talebi-Kalaleh](mailto:talebika@ualberta.ca)

        Welcome to the **Moving Load Analysis** application using [OpenSees](http://opensees.berkeley.edu/).
        This app is designed to guide you through model creation, load application,
        analysis, and post-processing of results in a clear and systematic way.
        """
    )

    # -- Sidebar for Navigation
    st.sidebar.title("Navigation")
    pages = ["Introduction", "Model Setup", "Load Application", "Analysis", "Results"]
    selected_page = st.sidebar.radio("Select a section:", pages)

    # -- Page Routing
    if selected_page == "Introduction":
        introduction_page()

    elif selected_page == "Model Setup":
        model_setup_page()

    elif selected_page == "Load Application":
        load_application_page()

    elif selected_page == "Analysis":
        analysis_page()

    elif selected_page == "Results":
        results_page()


def introduction_page():
    st.header("Introduction")
    st.write(
        """
        In this section, we'll describe the background and objective of the moving load analysis.
        - **Purpose**: To demonstrate how traffic or other moving loads affect a structural system.
        - **Software**: Built on top of OpenSees for finite element analysis.
        - **Workflow**: Model definition → Loading → Analysis → Post-processing.
        """
    )

def model_setup_page():
    st.header("Model Setup")

    st.write("### 1. Define Model Parameters")
    # Example form to collect model parameters from the user
    with st.form("model_parameters_form"):
        span_length = st.number_input("Span Length (m)", value=20.0, step=0.1)
        num_elements = st.number_input("Number of Elements", value=10, step=1)
        submitted = st.form_submit_button("Set Model Parameters")

    if submitted:
        st.session_state["model_params"] = {
            "span_length": span_length,
            "num_elements": num_elements,
        }
        st.success("Model parameters updated successfully!")

    st.write("### 2. Initialize the OpenSees Model")
    if st.button("Initialize Model"):
        # Retrieve model parameters from session state
        if "model_params" not in st.session_state:
            st.warning("Please set the model parameters first!")
        else:
            st.write("Initializing the model with parameters:", st.session_state["model_params"])
            create_opensees_model(st.session_state["model_params"])
            st.success("OpenSees model initialized successfully!")


def load_application_page():
    st.header("Load Application")

    st.write(
        """
        In this step, we define and apply the moving loads. The load could be
        modeled as a time-varying point load or multiple moving loads
        sequentially across the structure.
        """
    )

    if st.button("Apply Moving Load"):
        if "model_params" not in st.session_state:
            st.warning("Please define and initialize the model first!")
        else:
            apply_moving_load(st.session_state["model_params"])
            st.success("Moving load has been applied!")


def analysis_page():
    st.header("Analysis")

    st.write(
        """
        With the model and loads defined, we can now run the analysis. This could
        be a dynamic or quasi-static stepping analysis, depending on the approach.
        """
    )

    if st.button("Run Analysis"):
        if "model_params" not in st.session_state:
            st.warning("Please define and initialize the model first!")
        else:
            run_analysis(st.session_state["model_params"])
            st.success("Analysis completed!")


def results_page():
    st.header("Results")

    st.write(
        """
        Post-process and visualize the results of the moving load analysis here.
        You can view displacements, reactions, or any other relevant response quantities.
        """
    )

    col1, col2 = st.columns([1, 1])
    with col1:
        if st.button("Fetch Post-Processing Results"):
            results_data = post_processing_results()
            st.session_state["results_data"] = results_data
            st.success("Results data retrieved!")

    with col2:
        if st.button("Plot Results"):
            if "results_data" in st.session_state:
                plot_results(st.session_state["results_data"])
                st.info("Results plotted. Check your figures or console output.")
            else:
                st.warning("Please run post-processing to get results data first!")

    # Optionally, display a table or numeric results here
    if "results_data" in st.session_state:
        st.subheader("Displacements")
        st.write(st.session_state["results_data"].get("displacements", []))

        st.subheader("Reactions")
        st.write(st.session_state["results_data"].get("reactions", []))


if __name__ == "__main__":
    main()
