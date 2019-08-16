
from pele.storage import Database
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
import pylab as plt
import sys

class orderByValue(object):
      def __call__(self, min1):
         print 'ID',min1.coords[0]
         return min1.coords[0]

args = sys.argv
db_amp = Database(db=args[1],accuracy=1e-10)
db_lj = Database(db=args[2], accuracy=1e-10)
dg_amp = DisconnectivityGraph(database2graph(db_amp), nlevels=100, Emax=-56.7716+0.05,
                          order_by_basin_size=False,
                          center_gmin=False,
                          order_by_value=orderByValue())
dg_lj = DisconnectivityGraph(database2graph(db_lj), nlevels=100, Emax=-56.7716+0.05,
                          order_by_basin_size=True,
                          center_gmin=True,
                          order_by_value=orderByValue())
dg_lj.calculate()

fig = plt.figure(figsize=(9,7))
fig.set_facecolor('white')
ax = fig.add_subplot(111, adjustable='box')
#ax = fig.add_subplot(111)
#ax.axvline(linewidth=4, color="g")
ax.spines['left'].set_linewidth(10.0)
dg_lj.plot(linewidth=2., axes=ax, ylinewidth=1.5, xmin_offset=0.15)
#ylabel='Potential Energy')

plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.1)

dg_lj.draw_minima(db_lj.minima(), marker='o',c='black', s=50.0, linewidths=2.0, alpha=1)

dg_amp.calculate()
colors=[]
ms=[]

for m in db_amp.minima()[0:]:
   colors.append('red')
   ms.append([m])
print len(colors)
dg_amp.color_by_group(ms, colors=colors)
dg_amp.plot(linewidth=1.5, axes=ax, ylinewidth=1.5, xmin_offset=0.15)

emphasize_states = [3, 15]
emphasize = []
common_minima = []
for m in db_amp.minima():
  if int(m.coords[0]) in emphasize_states:
    emphasize.append(m)
  else:
    common_minima.append(m)
dg_amp.draw_minima(emphasize, marker='8',c='tab:green', s=50.0)
dg_amp.draw_minima(common_minima, marker='v',c='tab:red', s=20.0)

plt.yticks(fontsize=20.0, weight=40)
plt.tick_params(width=3, length=6.0)
plt.ylabel('Potential Energy', fontsize=22.0, fontweight='bold')
plt.savefig('./n7.png')
plt.show()
