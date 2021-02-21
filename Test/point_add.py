import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.widgets import SpanSelector


fig, ax = plt.subplots(figsize=(8, 6))
fig.subplots_adjust(bottom=0.2)

t = np.arange(0.0, 5.0, 0.01)
y = np.sin(2*np.pi*t) + 0.5*np.random.randn(len(t))

l = ax.plot(t, y, lw=2)
ax.set_title('Press left mouse button and drag to test')

new_data = []
indices = []

def onselect(tmin, tmax):
    indmin, indmax = np.searchsorted(t, (tmin, tmax))
    indmax = min(len(t) - 1, indmax)

    ab = [indmin, indmax]
    indices.append(ab)
    # save
    #np.savetxt("text.dat", np.c_[thisx, thisy])

def submit(expression):
    """
    Update the plotted function to the new math *expression*.

    *expression* is a string using "t" as its independent variable, e.g.
    "t ** 3".
    """
    ydata = expression
    new_data.append(ydata)
    text_box.set_val("")

# set useblit True on gtkagg for enhanced performance
span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'))


axbox = fig.add_axes([0.2, 0.05, 0.5, 0.075])
text_box = TextBox(axbox, "Enter the corresponding\n wavelength here (in Angstrom)")
text_box.on_submit(submit)
#text_box.set_val("")  # Trigger `submit` with the initial string.
text_box.stop_typing()

axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)

plt.show()

print(new_data)
print(indices)