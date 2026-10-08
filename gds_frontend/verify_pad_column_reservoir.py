"""Streaming Pad pruning must equal sorting the complete route library."""
from random import Random

from pad_router import _PadColumnReservoir,_diverse_pad_columns


def main():
    random=Random(1603)
    for size in (0,1,20,3000):
        for short,extra in ((1,0),(12,0),(12,6),(24,12)):
            choices=[]
            reservoir=_PadColumnReservoir(short,extra)
            for index in range(size):
                column={'id':index,
                        'source_node':random.randrange(1,max(2,size//4)),
                        'complete_length_um':random.randrange(1,100),
                        'depth_um':random.randrange(1,100)}
                choices.append(column)
                reservoir.add(column)
            expected=[c['id'] for c in
                      _diverse_pad_columns(choices,short,extra)]
            actual=[c['id'] for c in reservoir.selected()]
            assert actual==expected,(size,short,extra,expected,actual)
    print('streaming Pad candidate reservoir exactly matches full sorting',
          flush=True)


if __name__=='__main__':
    main()
