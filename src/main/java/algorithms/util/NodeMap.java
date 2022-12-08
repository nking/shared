/**
// This is an extension of and edits to TLongLongHashMap from the trove4j library
// downloaded from https://bitbucket.org/trove4j/trove/overview
// 
//   extended to add arrays to hold more values for the same key, specifically
     for nodes having parent, left, right, value, color and size.
// 
// The original copyright is:
// Copyright (c) 2001, Eric D. Friedman All Rights Reserved.
// Copyright (c) 2009, Rob Eden All Rights Reserved.
// Copyright (c) 2009, Jeff Randall All Rights Reserved.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// @author Eric D. Friedman
// @author Rob Eden
// @author Jeff Randall
// @version $Id: _K__V_HashMap.template,v 1.1.2.16 2010/03/02 04:09:50 robeden Exp $
*/
package algorithms.util;

import gnu.trove.map.TLongLongMap;
import gnu.trove.set.*;
import gnu.trove.iterator.*;
import gnu.trove.impl.hash.*;
import gnu.trove.impl.HashFunctions;
import gnu.trove.set.hash.TLongHashSet;

import java.io.*;
import java.util.*;

/**
 * An open addressed Map implementation for long keys and 
 *    long, long, long, int, int, int values.
 * 
 * Note that the structure was built specifically to hold fields for a tree node,
 * so there are special rules used to help define that a value has not been
 * set.
 * 
 * The key is the unique key for an entry.
 *    parent is LONG.MIN_VALUE when not used.
 *    left is LONG.MIN_VALUE when not used.
 *    right is LONG.MIN_VALUE when not used.
 *    value, as soon as the entry exists, always has a used value.
 *    color, as soon as the entry exists, always has a used value.
 *    size, as soon as the entry exists, always has a used value.
 *
 * Note that equals and hashcode use key and value only.
 * 
 */
public class NodeMap extends TLongIntHash implements 
    Externalizable {
    
    static final long serialVersionUID = 12325L;

    /** the values of the map */
    protected transient long[] _parents;

    /**
     *
     */
    protected transient long[] _lefts;

    /**
     *
     */
    protected transient long[] _rights;

    /**
     *
     */
    protected transient int[] _values;

    /**
     *
     */
    protected transient int[] _colors;

    /**
     *
     */
    protected transient int[] _sizes;
    
    /**
     *
     */
    protected final long noLinkValue = Long.MIN_VALUE;

    /**
     * Creates a new <code>TLongLongHashMap</code> instance with the default
     * capacity and load factor.
     */
    public NodeMap() {
        super();
        no_entry_key = 0;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance with a prime
     * capacity equal to or greater than <em>initialCapacity</em> and
     * with the default load factor.
     *
     @param initialCapacity an <code>int</code> value
     */
    public NodeMap( int initialCapacity ) {
        super( initialCapacity );
        no_entry_key = 0;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance with a prime
     * capacity equal to or greater than <em>initialCapacity</em> and
     * with the specified load factor.
     *
     @param initialCapacity an <code>int</code> value
     @param loadFactor a <code>float</code> value
     */
    public NodeMap( int initialCapacity, float loadFactor ) {
        super( initialCapacity, loadFactor );
        no_entry_key = 0;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance with a prime
     * capacity equal to or greater than <em>initialCapacity</em> and
     * with the specified load factor.
     * 
     * Note that the "no link value" is fixed to Long.MIN_VALUE.
     *
     @param initialCapacity an <code>int</code> value
     @param loadFactor a <code>float</code> value
     @param noEntryKey a <code>long</code> value that represents
     *                   <em>null</em> for the Key set.
     @param noEntryValue a <code>int</code> value that represents
     *                   <em>null</em> for the Value set.
     */
    public NodeMap( int initialCapacity, float loadFactor,
        long noEntryKey, int noEntryValue) {
        super( initialCapacity, loadFactor, noEntryKey, noEntryValue);
        no_entry_key = 0;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance containing
     * all of the entries in the map passed in.
     *      *
     @param keys a <em>long</em> array containing the keys for the matching values.
     *  
     @param parents is <em>long</em> array containing the first entry of values.
     *     <em>For trees using this class</em>, this is the <em>parent</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param lefts is <code>long</code> array containing the second entry of values.
     *     <em>For trees using this class</em>, this is the <em>left</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param rights is <code>long</code> array containing the third entry of values.
     *     <em>For trees using this class</em>, this is the <em>right</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param values is <code>int</code> array containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param colors is <code>int</code> array containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param sizes <code>int</code> array containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
     */
    public NodeMap( long[] keys, 
        long[] parents, long[] lefts, long[] rights, 
        int[] values, int[] colors, int[] sizes) {
        
        super(Math.max(keys.length, parents.length ) );
        
        no_entry_key = 0;
        
        if (keys.length != parents.length ||
            keys.length != lefts.length ||
            keys.length != rights.length ||
            keys.length != values.length ||
            keys.length != colors.length ||
            keys.length != sizes.length) {
            throw new IllegalArgumentException("all arrays must be same length");
        }

        int size = Math.min(keys.length, parents.length );
        for ( int i = 0; i < size; i++ ) {
            this.put(keys[i], 
                parents[i], lefts[i], rights[i],
                values[i], colors[i], sizes[i]);
        }
    }
    
    /**
     * from THash, a fast ceiling 
     @param v
     @return 
     */
    protected static long fastCeil( double v ) {
        long possible_result = ( long ) v;
        if ( v - possible_result > 0 ) possible_result++;
        return possible_result;
    }

    /* precondition: v > 0 */

    /**
     *
     @param v
     @return
     */

    protected static int saturatedCast(long v) {
        int r = (int) (v & 0x7fffffff); // removing sign bit
        if (r != v) {
            return Integer.MAX_VALUE;
        }
        return r;
    }

    /**
     * Creates a new <code>TLongLongHashMap</code> instance containing
     * all of the entries in the map passed in.
     *
     @param map a <code>TLongLongMap</code> that will be duplicated.
     */
    public NodeMap( 
        NodeMap map ) {
        
        super( map.size() );
        
        NodeMap hashmap 
            = ( NodeMap ) map;
        this._loadFactor = Math.abs( hashmap._loadFactor );
        this.no_entry_key = hashmap.no_entry_key;
        this.no_entry_value = hashmap.no_entry_value;
        
        //noinspection RedundantCast
        if ( this.no_entry_key != ( long ) 0 ) {
            Arrays.fill( _set, this.no_entry_key );
        }
                
        setUp( saturatedCast( 
            fastCeil( DEFAULT_CAPACITY / (double) _loadFactor ) ) );
         
        Arrays.fill(_parents, this.noLinkValue );
        Arrays.fill(_lefts, this.noLinkValue );
        Arrays.fill(_rights, this.noLinkValue );
        Arrays.fill(_values, this.no_entry_value );
        Arrays.fill(_colors, this.no_entry_value );
        Arrays.fill(_sizes, this.no_entry_value );
       
        putAll( map );
    }


    /**
     * initializes the hashtable to a prime capacity which is at least
     * <code>initialCapacity + 1</code>.
     *
     @param initialCapacity an <code>int</code> value
     @return the actual capacity chosen
     */
    protected int setUp( int initialCapacity ) {
        int capacity;

        capacity = super.setUp( initialCapacity );
        _parents = new long[capacity];
        _lefts = new long[capacity];
        _rights = new long[capacity];
        _values = new int[capacity];
        _colors = new int[capacity];
        _sizes = new int[capacity];
        return capacity;
    }


    /**
     * rehashes the map to the new capacity.
     *
     @param newCapacity an <code>int</code> value
     */
     /** {@inheritDoc} */
    protected void rehash( int newCapacity ) {
        
        int oldCapacity = _set.length;
        
        long oldKeys[] = _set;
        long oldVals0[] = _parents;
        long oldVals1[] = _lefts;
        long oldVals2[] = _rights;
        int oldVals3[] = _values;
        int oldVals4[] = _colors;
        int oldVals5[] = _sizes;
        byte oldStates[] = _states;

        _set = new long[newCapacity];
        _parents = new long[newCapacity];
        _lefts = new long[newCapacity];
        _rights = new long[newCapacity];
        _values = new int[newCapacity];
        _colors = new int[newCapacity];
        _sizes = new int[newCapacity];
        _states = new byte[newCapacity];

        for ( int i = oldCapacity; i-- > 0; ) {
            if( oldStates[i] == FULL ) {
                long o = oldKeys[i];
                int index = insertKey( o );
                _parents[index] = oldVals0[i];
                _lefts[index] = oldVals1[i];
                _rights[index] = oldVals2[i];
                _values[index] = oldVals3[i];
                _colors[index] = oldVals4[i];
                _sizes[index] = oldVals5[i];
            }
        }
    }

    /**
     * insert key and values.
     * 
     @param key a <code>long</code> containing the key for the associated values.
     *  
     @param parent is <code>long</code> containing the first entry of values.
     *     <em>For trees using this class</em>, this is the <em>parent</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param left is <code>long</code> containing the second entry of values.
     *     <em>For trees using this class</em>, this is the <em>left</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param right is <code>long</code> containing the third entry of values.
     *     <em>For trees using this class</em>, this is the <em>right</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param nodeValue is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeColor is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeSize <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
     */
    public void put( long key, long parent, long left, long right,
        int nodeValue, int nodeColor, int nodeSize) {
        int index = insertKey( key );
        doPut(key, parent, left, right, nodeValue, nodeColor, nodeSize, index );
    }
    
    /**
     insert key and values.  Note that the missing value0, value1, and value2
     are inserted with the sentinel for long null values, Long.MIN_VALUE.
     * 
     @param key a <code>long</code> containing the key for the associated values.
     *  
       @param nodeValue is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeColor is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeSize <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields. 
     */
    public void put( long key, int nodeValue, int nodeColor, int nodeSize) {
        int index = insertKey( key );
        doPut(key, noLinkValue, noLinkValue, noLinkValue, nodeValue, nodeColor, 
            nodeSize, index );
    }

    /** insert into map if there is not already an existing key.
     * Note that if there is an existing key, false is returned,
     * else true for successful insert.
     @param key a <code>long</code> containing the key for the associated values.
     *  
     @param parent is <code>long</code> containing the first entry of values.
     *     <em>For trees using this class</em>, this is the <em>parent</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param left is <code>long</code> containing the second entry of values.
     *     <em>For trees using this class</em>, this is the <em>left</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param right is <code>long</code> containing the third entry of values.
     *     <em>For trees using this class</em>, this is the <em>right</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param nodeValue is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeColor is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeSize <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
     @return 
     */
    public boolean putIfAbsent( long key, long parent, long left, long right,
        int nodeValue, int nodeColor, int nodeSize ) {
        int index = insertKey( key );
        if (index < 0)
            return false;
        doPut(key, parent, left, right, nodeValue, nodeColor, nodeSize, index );
        return true;
    }
    
    /** insert into map if there is not already an existing key.
     * Note that if there is an existing key, false is returned,
     * else true for successful insert.
     * 
     * Note that the missing value0, value1, and value2
     are inserted with the sentinel for long null values, Long.MIN_VALUE.
     * 
     @param key a <code>long</code> containing the key for the associated values.
     *  
       @param nodeValue is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeColor is <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param nodeSize <code>int</code> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields. 
     @return  
     */
    public boolean putIfAbsent( long key, int nodeValue, int nodeColor, int nodeSize ) {
        int index = insertKey( key );
        if (index < 0)
            return false;
        doPut(key, noLinkValue, noLinkValue, noLinkValue, nodeValue, nodeColor, nodeSize, index );
        return true;
    }

    private void doPut( long key, long parent, long left, long right,
        int value, int nodeColor, int nodeSize, int index ) {
        
        boolean isNewMapping = true;
        if ( index < 0 ) {
            index = -index -1;
            isNewMapping = false;
        }
        _parents[index] = parent;
        _lefts[index] = left;
        _rights[index] = right;
        _values[index] = value;
        _colors[index] = nodeColor;
        _sizes[index] = nodeSize;

        if (isNewMapping) {
            // a method in THash to expand and rehash if needed
            postInsertHook( consumeFreeSlot );
        }
    }
    
    /**
     * unset parent if it exists, else return false because there is
     * not an entry already for the given key.
     @param key
     @return 
     */
    public boolean unsetParent( long key) {
        //for existing entries, index returned by insertKey is negative
        int index = insertKey( key );
        if (index >= 0)
            return false;
        if ( index < 0 ) {
            index = -index -1;
        }
        _parents[index] = noLinkValue;
        return true;
    }
    /**
     * unset left if it exists, else return false because there is
     * not an entry already for the given key.
     @param key
     @return 
     */
    public boolean unsetLeft( long key) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _lefts[index] = noLinkValue;
        return true;
    }
    /**
     * unset right if it exists, else return false because there is
     * not an entry already for the given key.
     @param key
     @return 
     */
    public boolean unsetRight( long key ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _rights[index] = noLinkValue;
        return true;
    }
    
    /**
     * update parent if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     @param key
     @param parent
     @return 
     */
    public boolean updateParent( long key, long parent ) {
        //for existing entries, index returned by insertKey is negative
        int index = insertKey( key );
        if (index >= 0)
            return false;
        if ( index < 0 ) {
            index = -index -1;
        }
        _parents[index] = parent;
        return true;
    }
    /**
     * update left if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     @param key
     @param left
     @return 
     */
    public boolean updateLeft( long key, long left ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _lefts[index] = left;
        return true;
    }
    /**
     * update right if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     @param key
     @param right
     @return 
     */
    public boolean updateRight( long key, long right ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _rights[index] = right;
        return true;
    }
    /**
     * update nodeValue if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     @param key
     @param nodeValue
     @return 
     */
    public boolean updateNodeValue( long key, int nodeValue ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _values[index] = nodeValue;
        return true;
    }
    /**
     * update nodeColor if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     @param key
     @param nodeColor
     @return 
     */
    public boolean updateNodeColor( long key, int nodeColor ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _colors[index] = nodeColor;
        return true;
    }
    /**
     * update nodeSize if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     @param key
     @param nodeSize the size of the subtree underneath and including the node
     @return 
     */
    public boolean updateNodeSize( long key, int nodeSize ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _sizes[index] = nodeSize;
        return true;
    }
    /** {@inheritDoc} */
    public void putAll( NodeMap map ) {
        ensureCapacity( map.size() );
        TLong_LongLongLongIntIntIntHashIterator iter 
            = map.iterator();
        while ( iter.hasNext() ) {
            iter.advance();
            this.put( iter.key(), 
                iter.parent(), iter.left(), iter.right(),
                iter.nodeValue(), iter.nodeColor(), iter.nodeSize());
        }
    }

    /**
     *
     @param key
     @return
     */
    public boolean parentIsSet(long key) {
        int index = index(key);
        if (index < 0 || _parents[index] == noLinkValue) {
            return false;
        }
        return true;
    }

    /**
     *
     @param key
     @return
     */
    public boolean leftIsSet(long key) {
        int index = index(key);
        if (index < 0 || _lefts[index] == noLinkValue) {
            return false;
        }
        return true;
    }

    /**
     *
     @param key
     @return
     */
    public boolean rightIsSet(long key) {
        int index = index(key);
        if (index < 0 || _rights[index] == noLinkValue) {
            return false;
        }
        return true;
    }
    /** {@inheritDoc} */
    public long getParent( long key ) {
        int index = index( key );
        return index < 0 ? noLinkValue : _parents[index];
    }

    /**
     *
     @param key
     @return
     */
    public long getLeft( long key ) {
        int index = index( key );
        return index < 0 ? noLinkValue : _lefts[index];
    }

    /**
     *
     @param key
     @return
     */
    public long getRight( long key ) {
        int index = index( key );
        return index < 0 ? noLinkValue : _rights[index];
    }

    /**
     *
     @param key
     @return
     */
    public int getNodeValue( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value : _values[index];
    }

    /**
     *
     @param key
     @return
     */
    public int getNodeColor( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value : _colors[index];
    }

    /**
     *
     @param key
     @return
     */
    public int getNodeSize( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value : _sizes[index];
    }


    /** {@inheritDoc} */
    @Override
    public void clear() {
        super.clear();
        int n = _parents.length;
        Arrays.fill( _set, 0, _set.length, no_entry_key );
        Arrays.fill(_parents, 0, n, noLinkValue );
        Arrays.fill(_lefts, 0, n, noLinkValue );
        Arrays.fill(_rights, 0, n, noLinkValue );
        Arrays.fill(_values, 0, n, no_entry_value );
        Arrays.fill(_colors, 0, n, no_entry_value );
        Arrays.fill(_sizes, 0, n, no_entry_value );
        Arrays.fill( _states, 0, _states.length, FREE );
    }


    /** {@inheritDoc} */
    public boolean isEmpty() {
        return 0 == _size;
    }


    /** remove the entry for key and return true if succeeded, else false if
     * entry with key was not in map.
     @param key
     @return 
     */
    public boolean remove( long key ) {
        int index = index( key );
        if ( index >= 0 ) {
            removeAt( index );    // clear key,state; adjust size
            return true;
        }
        return false;
    }


    /** {@inheritDoc} */
    protected void removeAt( int index ) {
        _parents[index] = noLinkValue;
        _lefts[index] = noLinkValue;
        _rights[index] = noLinkValue;
        _values[index] = no_entry_value;
        _colors[index] = no_entry_value;
        _sizes[index] = no_entry_value;
        super.removeAt( index );  // clear key, state; adjust size
    }


    /** {@inheritDoc} */
    public TLongSet keySet() {
        
        TLongSet keys = new TLongHashSet();
        
        if (size() == 0) {
            return keys;
        }
        
        long[] k = _set;
        byte[] states = _states;

        for (int i = k.length, j = 0; i-- > 0;) {
            if (states[i] == FULL) {
                keys.add(k[i]);
            }
        }
        
        return keys;
    }


    /** {@inheritDoc} */
    public long[] keys() {
        long[] keys = new long[size()];
        if ( keys.length == 0 ) {
            return keys;        // nothing to copy
        }
        long[] k = _set;
        byte[] states = _states;

        for ( int i = k.length, j = 0; i-- > 0; ) {
          if ( states[i] == FULL ) {
            keys[j++] = k[i];
          }
        }
        return keys;
    }


    /** {@inheritDoc} */
    public long[] keys( long[] array ) {
        int size = size();
        if ( size == 0 ) {
            return array;       // nothing to copy
        }
        if ( array.length < size ) {
            array = new long[size];
        }

        long[] keys = _set;
        byte[] states = _states;

        for ( int i = keys.length, j = 0; i-- > 0; ) {
          if ( states[i] == FULL ) {
            array[j++] = keys[i];
          }
        }
        return array;
    }

    /**
     * same as contains key
     @param val
     @return 
     */
    @Override
    public boolean contains(long val) {
        return super.contains(val); 
    }
    
    /** {@inheritDoc} */
    public boolean containsKey( long key ) {
        return contains( key );
    }

    /**
     *
     @return
     */
    public TLong_LongLongLongIntIntIntHashIterator iterator() {
        return new TLong_LongLongLongIntIntIntHashIterator(this );
    }

    class TLong_LongLongLongIntIntIntHashKeyIterator 
        extends THashPrimitiveIterator implements TLongIterator {

        /**
         * Creates an iterator over the specified map
         *
         @param hash the <code>TPrimitiveHash</code> we will be iterating over.
         */
        TLong_LongLongLongIntIntIntHashKeyIterator( TPrimitiveHash hash ) {
            super( hash );
        }

        /** {@inheritDoc} */
        public long next() {
            moveToNextIndex();
            return _set[_index];
        }

        @Override
        public void remove() {
            if ( _expectedSize != _hash.size() ) {
                throw new ConcurrentModificationException();
            }

            // Disable auto compaction during the remove. This is a workaround for bug 1642768.
            try {
                _hash.tempDisableAutoCompaction();
                NodeMap.this.removeAt( _index );
            }
            finally {
                _hash.reenableAutoCompaction( false );
            }

            _expectedSize--;
        }
    }

    class TLong_LongLongLongIntIntIntHashIterator 
        extends THashPrimitiveIterator {

        /**
         * Creates an iterator over the specified map
         *
         @param map the <code>TLongLongHashMap</code> we will be iterating over.
         */
        TLong_LongLongLongIntIntIntHashIterator( 
            NodeMap map ) {
            super( map );
        }

        /** {@inheritDoc} */
        public void advance() {
            moveToNextIndex();
        }

        /** {@inheritDoc} */
        public long key() {
            return _set[_index];
        }

        /** {@inheritDoc} */
        public long parent() {
            return _parents[_index];
        }
        public long left() {
            return _lefts[_index];
        }
        public long right() {
            return _rights[_index];
        }
        public int nodeValue() {
            return _values[_index];
        }
        public int nodeColor() {
            return _colors[_index];
        }
        public int nodeSize() {
            return _sizes[_index];
        }

        /** {@inheritDoc} */
        public long setParent( long val ) {
            long old = parent();
            _parents[_index] = val;
            return old;
        }
        public long setLeft( long val ) {
            long old = left();
            _lefts[_index] = val;
            return old;
        }
        public long setRight( long val ) {
            long old = right();
            _rights[_index] = val;
            return old;
        }
        
        public int setValue( int val ) {
            int old = nodeValue();
            _values[_index] = val;
            return old;
        }
        public int setNodeColor( int val ) {
            int old = nodeColor();
            _colors[_index] = val;
            return old;
        }
        public int setNodeSize( int val ) {
            int old = nodeSize();
            _sizes[_index] = val;
            return old;
        }

        @Override
        public void remove() {
            if ( _expectedSize != _hash.size() ) {
                throw new ConcurrentModificationException();
            }
            // Disable auto compaction during the remove. This is a workaround for bug 1642768.
            try {
                _hash.tempDisableAutoCompaction();
                NodeMap.this.removeAt( _index );
            }
            finally {
                _hash.reenableAutoCompaction( false );
            }
            _expectedSize--;
        }
    }


    /** {@inheritDoc}
     * NOTE: this uses just the key and value3 for equals operation.
     */
    @Override
    public boolean equals( Object other ) {
        if ( ! ( other instanceof NodeMap ) ) {
            return false;
        }
        NodeMap that = ( NodeMap ) other;
        
        //System.out.println("in equals.  sizes=" + 
        //    that.size() + ", " + this.size());
        
        if ( that.size() != this.size() ) {
            return false;
        }
        int[] values3 = _values;
        byte[] states = _states;
      
        long this_no_entry_value = getNoEntryValue();
        long that_no_entry_value = that.getNoEntryValue();
        
        //System.out.println("this_no_entry_value=" + this_no_entry_value);
        //System.out.println("that_no_entry_value=" + that_no_entry_value);
        
        for ( int i = values3.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                long key = _set[i];

                //System.out.println("that.containsKey( key )="
                //    + that.containsKey(key));
                
                if ( !that.containsKey( key ) ) {
                    return false;
                }

                long that_value = that.getNodeValue(key);
                long this_value = values3[i];
                
                //System.out.println("values=" + that_value + ", " + this_value);
                
                if (
                    (this_value != that_value) && ( (this_value != this_no_entry_value)
                    || (that_value != that_no_entry_value))
                    ) {

                    return false;
                }
            }
        }
        return true;
    }


    /** {@inheritDoc} 
     * NOTE: this uses the key and values3 for hash code and no other fields.
     */
    @Override
    public int hashCode() {
        int hashcode = 0;
        byte[] states = _states;
        for ( int i = _parents.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                hashcode += HashFunctions.hash( _set[i] ) ^
                            HashFunctions.hash(_values[i] );
            }
        }
        return hashcode;
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {

        final StringBuilder buf = new StringBuilder( "{" );
        
        int[] values3 = _values;
        byte[] states = _states;
        for ( int i = values3.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                if (buf.length() > 0) {
                    buf.append( ", " );
                }
                long key = _set[i];
                buf.append(key);
                buf.append("=");
                buf.append(_parents[i]);
                buf.append(_lefts[i]);
                buf.append(_rights[i]);
                buf.append(_values[i]);
                buf.append(_colors[i]);
                buf.append(_sizes[i]);
            }
        }
        buf.append( "}" );
        return buf.toString();
    }

    /** {@inheritDoc} */
    public void writeExternal(ObjectOutput out) throws IOException {
        // VERSION
    	out.writeByte( 0 );

        // SUPER
    	super.writeExternal( out );

    	// NUMBER OF ENTRIES
    	out.writeInt( _size );

    	// ENTRIES
    	for ( int i = _states.length; i-- > 0; ) {
            if ( _states[i] == FULL ) {
                out.writeLong( _set[i] );
                out.writeLong(_parents[i] );
                out.writeLong(_lefts[i] );
                out.writeLong(_rights[i] );
                out.writeInt(_values[i] );
                out.writeInt(_colors[i] );
                out.writeInt(_sizes[i] );
            }
        }
        
        out.flush();
    
        //System.out.println("write map.  size=" + size());
    }


    /** {@inheritDoc} */
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
        
        // VERSION
    	in.readByte();

        // SUPER
    	super.readExternal( in );

    	// NUMBER OF ENTRIES
    	int size = in.readInt();
        
        //System.out.println("read size var=" + size);
    	
        setUp( size );
        
    	// ENTRIES
        for (int i = 0; i < size; ++i) {
            long key = in.readLong();
            long val0 = in.readLong();
            long val1 = in.readLong();
            long val2 = in.readLong();
            int val3 = in.readInt();
            int val4 = in.readInt();
            int val5 = in.readInt();
            
            put(key, val0, val1, val2, val3, val4, val5);
        }
        
        //System.out.println("read in map.  size=" + size);
    }
    
    /**
     *
     @param numberOfEntries
     @return
     */
    public static long estimateSizeOnHeap(int numberOfEntries) {
        
        /*
        super class:
           long serialVersionUID = 1L;
           long[] _set;
           long no_entry_key;
           int no_entry_value;
           boolean consumeFreeSlot;
         this class:
           long serialVersionUID = 12325L;
           long[] _parents;
           long[] _lefts;
           long[] _rights;
           int[] _values;
           int[] _colors;
           int[] _sizes;
           long noLinkValue = Long.MIN_VALUE;
        */
        
        ObjectSpaceEstimator est0 = new ObjectSpaceEstimator();
        est0.setNLongFields(2 + numberOfEntries);
        est0.setNIntFields(1);
        est0.setNArrayRefsFields(1);
        est0.setNBooleanFields(1);
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNLongFields(2 + 3*numberOfEntries);
        est.setNIntFields(1 + 3*numberOfEntries);
        est.setNArrayRefsFields(6);
        
        return est0.estimateSizeOnHeap() + est.estimateSizeOnHeap();
        
     }
} 
