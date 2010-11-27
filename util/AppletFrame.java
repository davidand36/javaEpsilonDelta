/*
  AppletFrame.java

  Wraps an applet in a JFrame, allowing a program to run as either an applet
  (on a Web page) or a stand-alone application.
  From Cay Horstmann & Gary Cornell, "Core Java", Vol I, Ch. 10.
*/

package us.EpsilonDelta.util;

import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.applet.*;
import java.io.*;
import java.net.*;


//*****************************************************************************


public
class AppletFrame
    extends JFrame
    implements AppletStub, AppletContext
{                                                                 //AppletFrame
//-----------------------------------------------------------------------------

    public
    AppletFrame( Applet applet )
    {
        m_applet = applet;
        add( m_applet );
        m_applet.setStub( this );
    }

//=============================================================================

    public
    void
    setVisible( boolean visible )
    {
        if ( visible )
        {
            m_applet.init( );
            super.setVisible( true );
            m_applet.start( );
        }
        else
        {
            m_applet.stop( );
            super.setVisible( false );
            m_applet.destroy( );
        }
    }

//=============================================================================

    //AppletStub methods

    public
    AppletContext
    getAppletContext( )
    {
        return this;
    }

//-----------------------------------------------------------------------------

    public
    boolean
    isActive( )
    {
        return true;
    }

//-----------------------------------------------------------------------------

    public
    URL
    getDocumentBase( )
    {
        return null;
    }

//-----------------------------------------------------------------------------

    public
    URL
    getCodeBase( )
    {
        return null;
    }

//-----------------------------------------------------------------------------

    public
    String
    getParameter( String name )
    {
        return "";
    }

//-----------------------------------------------------------------------------

    public
    void
    appletResize( int width, int height )
    {
    }

//=============================================================================

    //AppletContext methods

    public
    AudioClip
    getAudioClip( URL url )
    {
        return null;
    }

//-----------------------------------------------------------------------------

    public
    Image
    getImage( URL url )
    {
        return Toolkit.getDefaultToolkit().getImage( url );
    }

//-----------------------------------------------------------------------------

    public
    Applet
    getApplet( String name )
    {
        return null;
    }

//-----------------------------------------------------------------------------

    public
    Enumeration< Applet >
    getApplets( )
    {
        return null;
    }

//-----------------------------------------------------------------------------

    public
    void
    showDocument( URL url )
    {
    }

//.............................................................................

    public
    void
    showDocument( URL url, String target )
    {
    }
    
//-----------------------------------------------------------------------------

    public
    void
    showStatus( String status )
    {
    }

//-----------------------------------------------------------------------------

    public
    void
    setStream( String key, InputStream stream )
    {
    }

//-----------------------------------------------------------------------------

    public
    InputStream
    getStream( String key )
    {
        return null;
    }

//-----------------------------------------------------------------------------

    public
    Iterator< String >
    getStreamKeys( )
    {
        return null;
    }
    
//=============================================================================

    private Applet m_applet;
    
//-----------------------------------------------------------------------------
}                                                                 //AppletFrame


//*****************************************************************************
