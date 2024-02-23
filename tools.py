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
