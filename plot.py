import plotly.express as px


def plot_data(time, bar, sample):
    # BAR FORCE TRANSMITTED
    try:
        fig = px.line(x=time, y=bar['transmitted force'], title='Force transmitted')
        fig.show()
    except:
        pass

    # BAR FORCE REFLECTED
    try:
        fig = px.line(x=time, y=bar['reflected force'], title='Force reflected')
        fig.show()
    except:
        pass

    # SAMPLE STRESS
    try:
        fig = px.line(x=time, y=sample['stress']/1e6, title='Sample Stress')
        fig.show()
    except:
        pass

    # SAMPLE STRAIN
    try:
        fig = px.line(x=time, y=sample['strain'], title='Sample Strain')
        fig.show()
    except:
        pass

    # SAMPLE STRAIN RATE
    try:
        fig = px.line(x=time, y=sample['strain-rate'], title='Sample Strain Rate')
        fig.show()
    except:
        pass

    # STRESS-STRAIN ENGINEERING NOT CORRECTED
    try:
        fig = px.line(x=sample['strain'], y=sample['stress']/1e6, title='Stress-Strain Engineering not corrected')
        fig.show()
    except:
        pass

    # STRESS-STRAIN ENGINEERING CORRECTED
    try:
        fig = px.line(x=sample['strain extended'], y=sample['stress extended']/1e6, title='Stress-Strain Engineering corrected')
        fig.show()
    except:
        pass

    # STRESS-STRAIN TRUE
    try:
        fig = px.line(x=sample['strain true'], y=sample['stress true']/1e6, title='Stress-Strain True')
        fig.show()
    except:
        pass

