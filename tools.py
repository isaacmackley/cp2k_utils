def lighten_hex_colour(hex_colour, factor=0.2):
    # Convert hex colour to RGB components
    r = int(hex_colour[1:3], 16)
    g = int(hex_colour[3:5], 16)
    b = int(hex_colour[5:7], 16)
    
    # Increase the values of RGB components to lighten the colour
    r = min(255, int(r * (1 + factor)))
    g = min(255, int(g * (1 + factor)))
    b = min(255, int(b * (1 + factor)))
    
    # Convert the modified RGB values back to hexadecimal colour code
    return "#{:02X}{:02X}{:02X}".format(r, g, b)


def style_axes(ax, ticks_fontsize=20, hide_y=False):
    """
    House axis style shared by the CTL and PDOS plots: inward major/minor
    ticks, thicker spines, no grid. With hide_y=True the y tick marks and
    labels are removed (e.g. DOS in arbitrary units).
    """
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', direction='in',
                   length=6, width=1.5, labelsize=ticks_fontsize)
    ax.tick_params(axis='both', which='minor', direction='in', length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.grid(visible=False)
    if hide_y:
        ax.set_yticks([])
        ax.set_yticks([], minor=True)
