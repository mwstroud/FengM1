from bmtk.builder.networks import NetworkBuilder

net = NetworkBuilder('hco_net')
net.add_nodes(cell_name='HCOCell',
              model_type='biophysical',
              model_template='hoc:HCOcell',
              morphology='blank.swc',
              HCOCell='HCOCell'
            )
net.build()
net.save_nodes(output_dir='network')
