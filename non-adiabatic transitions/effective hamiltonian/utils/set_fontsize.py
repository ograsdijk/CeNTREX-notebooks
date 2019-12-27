def set_fontsize(ax, fs, legend = False):
    """
    Small function to set font sizes for figures
    """
    ax.tick_params(axis='both', which='major', labelsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)
    ax.xaxis.label.set_size(fs)
    ax.yaxis.label.set_size(fs)
    ax.yaxis.offsetText.set_fontsize(fs)
    if legend:
        [t.set_fontsize(fs) for t in ax.legend().get_texts()]
    ax.title.set_size(int(fs*1.2))