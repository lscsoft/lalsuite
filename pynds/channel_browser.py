#!/usr/bin/env python
"""

Channel browser for NDS1/NDS2 connections

"""
__author__       = "Leo Singer <leo.singer@ligo.org>"
__organization__ = ["LIGO", "California Institute of Technology"]
__copyright__    = "Copyright 2010, Leo Singer"




import gtk
import gobject
import nds


def split_channel(channel):
    channel_name = channel.name
    name_parts = []
    for i in range(len(channel_name)):
        if channel_name[i] in '-:_':
            name_parts.append(channel_name[:i])
        if len(name_parts) == 3:
            break
    name_parts.extend([name_parts[-1]] * (3 - len(name_parts)))
    name_parts.extend([channel.name, channel])
    return tuple(name_parts)


class ChannelTreeElement(object):
    def __init__(self, value=None):
        self.children = []
        self.selected = False
        self.value = value


class ChannelTreeNode(ChannelTreeElement):
    def __init__(self, value=None, parent=None):
        super(ChannelTreeNode, self).__init__(value)
        self.parent = parent
        self.countEnabledChildren = 0
        if parent is not None:
            self.parent.children.append(self)
    def incrementEnabledChildren(self):
        self.countEnabledChildren += 1
        if self.countEnabledChildren == 1:
            self.parent.incrementEnabledChildren()
    def decrementEnabledChildren(self):
        self.countEnabledChildren -= 1
        if self.countEnabledChildren == 0:
            self.parent.decrementEnabledChildren()
    def isEnabled(self):
        return self.countEnabledChildren > 0 and self.parent.selected


class ChannelTreeRoot(ChannelTreeElement):
    
    def __init__(self):
        super(ChannelTreeRoot, self).__init__()
        self.selected = True
    
    def incrementEnabledChildren(self):
        pass
    
    def decrementEnabledChildren(self):
        pass


class ChannelTreeLeaf(ChannelTreeNode):
    def incrementEnabledChildren(self):
        if self.countEnabledChildren < 1:
            super(ChannelTreeLeaf, self).incrementEnabledChildren()
    
    def decrementEnabledChildren(self):
        if self.countEnabledChildren > 0:
            super(ChannelTreeLeaf, self).decrementEnabledChildren()


def make_channel_tree(channels):
    tree = ChannelTreeRoot()
    instrument_node = ChannelTreeNode()
    leaves = list()
    
    for channel_parts in sorted(split_channel(c) for c in channels):
        instrument, system, subsystem, channel, channel_record = channel_parts
        
        if instrument != instrument_node.value:
            instrument_node = ChannelTreeNode(instrument, tree)
            system_node = ChannelTreeNode()
        
        if system != system_node.value:
            system_node = ChannelTreeNode(system, instrument_node)
            subsystem_node = ChannelTreeNode()
        
        if subsystem != subsystem_node.value:
            subsystem_node = ChannelTreeNode(subsystem, system_node)
        
        leaf = ChannelTreeLeaf(channel_record, subsystem_node)
        leaf.incrementEnabledChildren()
        leaves.append(leaf)
    
    return (tuple(leaves), tree)



class ChannelBrowser(gtk.HPaned):
    """An NDS1/NDS2 channel browser widget"""
    
    @staticmethod
    def list_store_add(liststore, objs):
        """Helper function to insert objects into a liststore"""
        for obj in objs:
            liststore.append( (obj,) )
    
    def __init__(self, channels):
        super(ChannelBrowser, self).__init__()
        self.channel_activated_func = None
        self.channelLeaves, self.channelTree = make_channel_tree(channels)
        self.rates = tuple(sorted(set(int(c.rate) for c in channels)))
        self.channel_types = tuple(c for c in nds.channel_type.values.values() if c != nds.channel_type.unknown)
        self.selected_rates = frozenset(self.rates)
        self.selected_channel_types = frozenset(self.channel_types)
        
        # Construct ListStore objects for four panes
        instrument_store = gtk.ListStore(gobject.TYPE_PYOBJECT)
        system_store     = gtk.ListStore(gobject.TYPE_PYOBJECT)
        subsystem_store  = gtk.ListStore(gobject.TYPE_PYOBJECT)
        channel_store    = gtk.ListStore(gobject.TYPE_PYOBJECT)
        
        # Populate the ListStore objects with elements from the channel tree
        self.list_store_add(instrument_store, self.channelTree.children)
        for instrument_node in self.channelTree.children:
            self.list_store_add(system_store, instrument_node.children)
            for system_node in instrument_node.children:
                self.list_store_add(subsystem_store, system_node.children)
                for subsystem_node in system_node.children:
                    self.list_store_add(channel_store, subsystem_node.children)
        
        # Construct the view
        scrolledwindow = gtk.ScrolledWindow()
        scrolledwindow.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        
        box1 = gtk.VBox()
        scrolledwindow.add_with_viewport(box1)
        self.pack1(scrolledwindow, shrink=False)
        
        liststore = gtk.ListStore(gobject.TYPE_UINT)
        self.list_store_add(liststore, self.rates)
        self.rates_filter = liststore.filter_new()
        self.rates_filter.set_modify_func((gobject.TYPE_STRING,), self.rates_filter_modify)
        treeview = gtk.TreeView(self.rates_filter)
        tvcolumn = gtk.TreeViewColumn('Rate (Hz)')
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText()
        tvcolumn.pack_start(cell, True)
        tvcolumn.add_attribute(cell, 'text', 0)
        treeview.set_search_column(0)
        treeview.set_reorderable(False)
        treeview.set_rubber_banding(True)
        selection = treeview.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)
        selection.select_all()
        selection.connect('changed', self.rates_selection_changed)
        box1.pack_start(treeview, expand=False, fill=False)
        
        liststore = gtk.ListStore(gobject.TYPE_PYOBJECT)
        self.list_store_add(liststore, self.channel_types)
        self.channel_types_filter = liststore.filter_new()
        self.channel_types_filter.set_modify_func((gobject.TYPE_STRING,), self.channel_types_filter_modify)
        treeview = gtk.TreeView(self.channel_types_filter)
        tvcolumn = gtk.TreeViewColumn('Channel type')
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText()
        tvcolumn.pack_start(cell, True)
        tvcolumn.add_attribute(cell, 'text', 0)
        treeview.set_search_column(0)
        treeview.set_reorderable(False)
        treeview.set_rubber_banding(True)
        selection = treeview.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)
        selection.select_all()
        selection.connect('changed', self.channel_types_selection_changed)
        box1.pack_start(treeview, expand=True, fill=True)
        
        box1 = gtk.HBox()
        box1.set_homogeneous(True)
        box1.set_size_request(640, 480)
        self.pack2(box1, resize=True, shrink=True)
        
        self.instruments_filter = instrument_store.filter_new()
        self.instruments_filter.set_visible_func(self.category_filter_visible)
        self.instruments_filter.set_modify_func((gobject.TYPE_STRING,), self.category_filter_modify)
        treeview = gtk.TreeView(self.instruments_filter)
        tvcolumn = gtk.TreeViewColumn('Instrument')
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText()
        tvcolumn.pack_start(cell, True)
        tvcolumn.add_attribute(cell, 'text', 0)
        treeview.set_search_column(0)
        treeview.set_reorderable(False)
        treeview.set_rubber_banding(True)
        selection = treeview.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)
        selection.set_select_function(self.category_select, full=True)
        selection.connect('changed', self.instruments_selection_changed)
        scrolledwindow = gtk.ScrolledWindow()
        scrolledwindow.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        scrolledwindow.add(treeview)
        box1.pack_start(scrolledwindow)
        
        self.systems_filter = system_store.filter_new()
        self.systems_filter.set_visible_func(self.category_filter_visible)
        self.systems_filter.set_modify_func((gobject.TYPE_STRING,), self.category_filter_modify)
        treeview = gtk.TreeView(self.systems_filter)
        tvcolumn = gtk.TreeViewColumn('System')
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText()
        tvcolumn.pack_start(cell, True)
        tvcolumn.add_attribute(cell, 'text', 0)
        treeview.set_search_column(0)
        treeview.set_reorderable(False)
        treeview.set_rubber_banding(True)
        selection = treeview.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)
        selection.set_select_function(self.category_select, full=True)
        selection.connect('changed', self.systems_selection_changed)
        scrolledwindow = gtk.ScrolledWindow()
        scrolledwindow.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        scrolledwindow.add(treeview)
        box1.pack_start(scrolledwindow)
        
        self.subsystems_filter = subsystem_store.filter_new()
        self.subsystems_filter.set_visible_func(self.category_filter_visible)
        self.subsystems_filter.set_modify_func((gobject.TYPE_STRING,), self.category_filter_modify)
        treeview = gtk.TreeView(self.subsystems_filter)
        tvcolumn = gtk.TreeViewColumn('Subsystem')
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText()
        tvcolumn.pack_start(cell, True)
        tvcolumn.add_attribute(cell, 'text', 0)
        treeview.set_search_column(0)
        treeview.set_reorderable(False)
        treeview.set_rubber_banding(True)
        selection = treeview.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)
        selection.set_select_function(self.category_select, full=True)
        selection.connect('changed', self.subsystems_selection_changed)
        scrolledwindow = gtk.ScrolledWindow()
        scrolledwindow.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        scrolledwindow.add(treeview)
        box1.pack_start(scrolledwindow)
        
        self.channels_filter = channel_store.filter_new()
        self.channels_filter.set_visible_func(self.category_filter_visible)
        self.channels_filter.set_modify_func((gobject.TYPE_STRING,), self.category_filter_modify)
        treeview = gtk.TreeView(self.channels_filter)
        tvcolumn = gtk.TreeViewColumn('Channel')
        treeview.connect('row-activated', self.channel_row_activated)
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText()
        tvcolumn.pack_start(cell, True)
        tvcolumn.add_attribute(cell, 'text', 0)
        treeview.set_search_column(0)
        treeview.set_reorderable(False)
        treeview.set_rubber_banding(True)
        selection = treeview.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)
        selection.set_select_function(self.category_select, full=True)
        scrolledwindow = gtk.ScrolledWindow()
        scrolledwindow.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
        scrolledwindow.add(treeview)
        box1.pack_start(scrolledwindow)
    
    def rates_filter_modify(self, model, iter, column):
        if column != 0:
            return None
        orig = model.get_model().get_value(model.convert_iter_to_child_iter(iter), 0)
        nMatches = 0
        for c in self.channelLeaves:
            channel = c.value
            if int(channel.rate) == orig and channel.type in self.selected_channel_types:
                nMatches += 1
        return '%d (%d)' % (orig, nMatches)
    
    def rates_selection_changed(self, object):
        self.selected_rates = frozenset(self.rates[p[0]] for p in object.get_selected_rows()[1])
        self.channel_types_filter.refilter()
        for channelLeaf in self.channelLeaves:
            if int(channelLeaf.value.rate) in self.selected_rates and channelLeaf.value.type in self.selected_channel_types:
                channelLeaf.incrementEnabledChildren()
            else:
                channelLeaf.decrementEnabledChildren()
        self.instruments_filter.refilter()
        self.systems_filter.refilter()
        self.subsystems_filter.refilter()
        self.channels_filter.refilter()
    
    def channel_types_filter_modify(self, model, iter, column):
        if column != 0:
            return None
        orig = model.get_model().get_value(model.convert_iter_to_child_iter(iter), 0)
        nMatches = 0
        for c in self.channelLeaves:
            channel = c.value
            if channel.type == orig and int(channel.rate) in self.selected_rates:
                nMatches += 1
        return '%s (%d)' % (str(orig), nMatches)
    
    def channel_types_selection_changed(self, object):
        self.selected_channel_types = frozenset(self.channel_types[p[0]] for p in object.get_selected_rows()[1])
        self.rates_filter.refilter()
        for channelLeaf in self.channelLeaves:
            if int(channelLeaf.value.rate) in self.selected_rates and channelLeaf.value.type in self.selected_channel_types:
                channelLeaf.incrementEnabledChildren()
            else:
                channelLeaf.decrementEnabledChildren()
        self.instruments_filter.refilter()
        self.systems_filter.refilter()
        self.subsystems_filter.refilter()
        self.channels_filter.refilter()
    
    @staticmethod
    def category_filter_visible(model, iter):
        return model.get_value(iter, 0).isEnabled()
    
    @staticmethod
    def category_filter_modify(model, iter, column):
        if column != 0:
            return None
        else:
            return str(model.get_model().get_value(model.convert_iter_to_child_iter(iter), 0).value)
    
    @staticmethod
    def category_select(selection, model, path, is_selected):
        model.get_model().get_value(model.convert_iter_to_child_iter(model.get_iter(path)), 0).selected = not(is_selected)
        return True
    
    def instruments_selection_changed(self, object):
        self.systems_filter.refilter()
        self.subsystems_filter.refilter()
        self.channels_filter.refilter()
    
    def systems_selection_changed(self, object):
        self.subsystems_filter.refilter()
        self.channels_filter.refilter()
    
    def subsystems_selection_changed(self, object):
        self.channels_filter.refilter()
    
    def channel_row_activated(self, treeview, path, view_column):
        if self.channel_activated_func:
            model = treeview.get_model()
            child_model = model.get_model()
            channel = child_model.get_value(model.convert_iter_to_child_iter(model.get_iter(path)), 0).value
            self.channel_activated_func(channel)



class Base:
    """A simple application that displays a channel browser in a window."""
    
    def __init__(self, host, port):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("destroy", self.destroy)
        
        daq = nds.daq(host, port)
        channels = daq.recv_channel_list()
        self.channelbrowser = ChannelBrowser(channels)
        self.window.set_title("%s:%d" % (daq.host, daq.port))
        
        self.channelbrowser.channel_activated_func = self.channel_activated
        self.window.add(self.channelbrowser)
        self.window.show_all()
    
    @staticmethod
    def channel_activated(channel):
        print channel
    
    @staticmethod
    def destroy(widget, data=None):
        gtk.main_quit()
    
    @staticmethod
    def main():
        gtk.main()



if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-n', '--host', default='blue.ligo-wa.caltech.edu', metavar='HOSTNAME')
    parser.add_option('-p', '--port', type='int', default=31200, metavar='PORT')
    (options, args) = parser.parse_args()
    
    base = Base(options.host, options.port)
    base.main()
