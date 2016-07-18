from array import array
import itertools

from0to13000 = (
array('f', [0.20333290100097656]),
array('f', [0.1911168247461319]),
array('f', [0.20421168208122253]),
array('f', [0.20437951385974884]),
array('f', [0.20393006503582]),
array('f', [0.0007833830895833671]),
array('f', [0.0007851485279388726]),
array('f', [0.0007813861593604088]),
array('f', [0.0007852490525692701]),
array('f', [0.0007824577041901648]),
array('f', [0.15376122295856476]),
array('f', [0.14818817377090454]),
array('f', [0.15722395479679108]),
array('f', [0.15372560918331146]),
array('f', [0.1543569713830948]),
array('f', [0.0008794961031526327]),
array('f', [0.0008835632470436394]),
array('f', [0.0008753624279052019]),
array('f', [0.0008824857068248093]),
array('f', [0.0008796491310931742]),
array('f', [0.3097667098045349]),
array('f', [0.30546635389328003]),
array('f', [0.31013211607933044]),
array('f', [0.25355425477027893]),
array('f', [0.29120442271232605]),
array('f', [0.0013285324675962329]),
array('f', [0.0013282749569043517]),
array('f', [0.0013287729816511273]),
array('f', [0.001328500802628696]),
array('f', [0.0013290469069033861]),
)

from105to140 = (
array('f', [0.20432309806346893]),
array('f', [0.19204752147197723]),
array('f', [0.20520614087581635]),
array('f', [0.2053748071193695]),
array('f', [0.20492316782474518]),
array('f', [0.032891370356082916]),
array('f', [0.03296549245715141]),
array('f', [0.032807525247335434]),
array('f', [0.03296971693634987]),
array('f', [0.03285251557826996]),
array('f', [0.15426567196846008]),
array('f', [0.14867433905601501]),
array('f', [0.1577397584915161]),
array('f', [0.15422993898391724]),
array('f', [0.1548633724451065]),
array('f', [0.031444206833839417]),
array('f', [0.03158961609005928]),
array('f', [0.031296417117118835]),
array('f', [0.031551092863082886]),
array('f', [0.03144967555999756]),
array('f', [0.312916100025177]),
array('f', [0.30857205390930176]),
array('f', [0.3132852613925934]),
array('f', [0.2561321556568146]),
array('f', [0.29416510462760925]),
array('f', [0.030636530369520187]),
array('f', [0.03063059225678444]),
array('f', [0.030642077326774597]),
array('f', [0.03063580021262169]),
array('f', [0.03064839355647564]),
)

wrong = {}
right = {}
for i, (channel, sorb, sys) in enumerate(itertools.product(("2e2mu", "4e", "4mu"), ("signal", "bkg"), ("", "ScaleUp", "ScaleDown", "ResUp", "ResDown"))):
    wrong[channel,sorb,sys] = from0to13000[i][0]
    right[channel,sorb,sys] = from105to140[i][0]

cconstantforpm4l = {}
for channel, sys in itertools.product(("2e2mu", "4e", "4mu"), ("", "ScaleUp", "ScaleDown", "ResUp", "ResDown")):
    cconstantforpm4l[channel, sys] = right[channel, "bkg", sys] / right[channel, "signal", sys] / (wrong[channel, "bkg", sys] / wrong[channel, "signal", sys])
    if __name__ == "__main__":
        print channel, sys, cconstantforpm4l[channel, sys]
    assert abs(cconstantforpm4l[channel, sys] / cconstantforpm4l[channel, ""] - 1) < 0.001

def cconstant(channel):
    if channel == 0: return cconstantforpm4l["4mu",""]
    if channel == 1: return cconstantforpm4l["4e",""]
    if channel == 2: return cconstantforpm4l["2e2mu",""]
    assert False
