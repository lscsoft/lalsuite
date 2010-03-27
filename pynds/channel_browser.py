#!/usr/bin/env python

import gtk
import gobject
import nds

def list_store_fill(liststore, iterable):
    for x in iterable:
        liststore.append((x,))
    return liststore


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



class Base:
    
    def __init__(self, host, port):
        self.daq = nds.daq(host, port)
        self.channels = self.daq.recv_channel_list()
        self.channelLeaves, self.channelTree = make_channel_tree(self.channels)
        self.rates = tuple(sorted(set(int(c.rate) for c in self.channels)))
        self.channel_types = tuple(c for c in nds.channel_type.values.values() if c != nds.channel_type.unknown)
        self.selected_rates = frozenset(self.rates)
        self.selected_channel_types = frozenset(self.channel_types)
    
    def rates_filter_modify(self, model, iter, column):
        if column != 0:
            return None
        orig = model.get_model().get_value(model.convert_iter_to_child_iter(iter), column)
        nMatches = 0
        for c in self.channels:
            if int(c.rate) == orig and c.type in self.selected_channel_types:
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
        orig = model.get_model().get_value(model.convert_iter_to_child_iter(iter), column)
        nMatches = 0
        for c in self.channels:
            if c.type == orig and int(c.rate) in self.selected_rates:
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
    
    def category_filter_visible(self, model, iter):
        return model.get_value(iter, 0).isEnabled()
    
    def category_filter_modify(self, model, iter, column):
        if column != 0:
            return None
        else:
            return str(model.get_model().get_value(model.convert_iter_to_child_iter(iter), 0).value)
    
    def category_select(self, selection, model, path, is_selected):
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
    
    def main(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("destroy", self.destroy)
        self.window.set_title("%s:%d" % (self.daq.host, self.daq.port))
        
        
        box = gtk.HPaned()
        
        scrolledwindow = gtk.ScrolledWindow()
        scrolledwindow.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        
        box1 = gtk.VBox()
        scrolledwindow.add_with_viewport(box1)
        box.pack1(scrolledwindow, shrink=False)
        
        liststore = list_store_fill(gtk.ListStore(gobject.TYPE_UINT), self.rates)
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
        
        liststore = list_store_fill(gtk.ListStore(gobject.TYPE_PYOBJECT), self.channel_types)
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
        box.pack2(box1, resize=True, shrink=True)
        
        liststore = list_store_fill(gtk.ListStore(gobject.TYPE_PYOBJECT), self.channelTree.children)
        self.instruments_filter = liststore.filter_new()
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
        
        liststore = gtk.ListStore(gobject.TYPE_PYOBJECT)
        for instrument_node in self.channelTree.children:
            list_store_fill(liststore, instrument_node.children)
        self.systems_filter = liststore.filter_new()
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
        
        liststore = gtk.ListStore(gobject.TYPE_PYOBJECT)
        for instrument_node in self.channelTree.children:
            for system_node in instrument_node.children:
                list_store_fill(liststore, system_node.children)
        self.subsystems_filter = liststore.filter_new()
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
        
        liststore = gtk.ListStore(gobject.TYPE_PYOBJECT)
        for instrument_node in self.channelTree.children:
            for system_node in instrument_node.children:
                for subsystem_node in system_node.children:
                    list_store_fill(liststore, subsystem_node.children)
        self.channels_filter = liststore.filter_new()
        self.channels_filter.set_visible_func(self.category_filter_visible)
        self.channels_filter.set_modify_func((gobject.TYPE_STRING,), self.category_filter_modify)
        treeview = gtk.TreeView(self.channels_filter)
        tvcolumn = gtk.TreeViewColumn('Channel')
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
        
        self.window.add(box)
        self.window.show_all()
        self.window.show()
        
        gtk.main()
    
    def destroy(self, widget, data=None):
        gtk.main_quit()




if __name__ == '__main__':
    base = Base("blue.ligo-wa.caltech.edu", 31200)
    base.main()
