
package main;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;
import javax.swing.JLabel;
import javax.swing.Timer;

public class ClockLabel extends JLabel implements ActionListener {
                                                            
        private SimpleDateFormat ft = new SimpleDateFormat ("E MMM.dd.yyyy - hh:mm:ss a", Locale.US);  
        private SimpleDateFormat ft2 = new SimpleDateFormat ("MMM-dd-yyyy (hh-mm-ss a)", Locale.US);
        
        private String alternateDate;
    
        public ClockLabel(){
            super("");
            super.setText(ft.format(new Date()));
            alternateDate = ft2.format(new Date());
            Timer timer = new Timer(1000, this);
            timer.start();
        }
        
        @Override
        public void actionPerformed(ActionEvent ae) {
            Date myDate = new Date( );
            setText(ft.format(myDate));
            alternateDate = ft2.format(new Date());
        }
        
        public String getDate(){
            return alternateDate;
        }
    }    